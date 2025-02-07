package id

import (
	"errors"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"path/filepath"
	"strings"

	"philosopher/lib/iso"
	"philosopher/lib/met"
	"philosopher/lib/msg"

	"philosopher/lib/mod"
	"philosopher/lib/spc"
	"philosopher/lib/sys"

	"github.com/sirupsen/logrus"
	"github.com/vmihailenco/msgpack"
)

// ProtXML struct
type ProtXML struct {
	FileName   string
	DecoyTag   string
	Groups     GroupList
	RunOptions string
}

// GroupIdentification tag
type GroupIdentification struct {
	GroupNumber uint32
	Probability float64
	Proteins    ProtIDList
}

// ProteinIdentification struct
type ProteinIdentification struct {
	GroupNumber              uint32
	GroupSiblingID           string
	ProteinName              string
	Description              string
	UniqueStrippedPeptides   []string
	Length                   string
	PercentCoverage          float32
	PctSpectrumIDs           float32
	GroupProbability         float64
	Probability              float64
	Confidence               float64
	TopPepProb               float64
	IndistinguishableProtein []string
	TotalNumberPeptides      int
	PeptideIons              []PeptideIonIdentification
	HasRazor                 bool
	Picked                   int
}

// PeptideIonIdentification struct
type PeptideIonIdentification struct {
	PeptideSequence          string
	ModifiedPeptide          string
	Charge                   uint8
	InitialProbability       float64
	Weight                   float64
	GroupWeight              float64
	CalcNeutralPepMass       float64
	NumberOfEnzymaticTermini uint8
	NumberOfInstances        int
	SharedParentProteins     int
	Razor                    int
	IsNondegenerateEvidence  bool
	IsUnique                 bool
	PeptideParentProtein     []string
	Labels                   iso.Labels
	Modifications            mod.Modifications
}

// GroupList represents a protein group list
type GroupList []GroupIdentification

// ProtIDList list represents a list of custom protein identifications
type ProtIDList []ProteinIdentification

// Len function for sortng
func (p ProtIDList) Len() int {
	return len(p)
}

// Less function for sorting
func (p ProtIDList) Less(i, j int) bool {
	return p[i].TopPepProb > p[j].TopPepProb
}

// Swap function for sorting
func (p ProtIDList) Swap(i, j int) {
	p[i], p[j] = p[j], p[i]
}

// Read ...
func (p *ProtXML) Read(f string) {

	var xml spc.ProtXML
	xml.Parse(f)

	var ptg = xml.ProteinSummary.ProteinGroup
	var groups GroupList

	for _, i := range ptg {

		var gi GroupIdentification
		var proteinList ProtIDList

		gi.GroupNumber = i.GroupNumber
		gi.Probability = i.Probability

		for jindex, j := range i.Protein {

			// correcting group probabilities
			if jindex == 0 {
				if i.Probability == 1 && j.Probability == 0 {
					j.Probability = float64(i.Probability)
					break
				}
			}

			var ptid ProteinIdentification

			ptid.GroupNumber = i.GroupNumber
			ptid.GroupProbability = i.Probability
			ptid.Probability = i.Probability
			ptid.ProteinName = string(j.ProteinName)
			ptid.Description = string(j.Annotation.ProteinDescription)
			ptid.Probability = j.Probability
			ptid.PercentCoverage = j.PercentCoverage
			ptid.PctSpectrumIDs = j.PctSpectrumIDs
			ptid.GroupSiblingID = string(j.GroupSiblingID)
			ptid.TotalNumberPeptides = j.TotalNumberPeptides
			ptid.TopPepProb = 0

			if strings.EqualFold(j.Parameter.Name, "prot_length") {
				ptid.Length = j.Parameter.Value
			}

			// collect indistinguishable proteins (Protein to Protein equivalency)
			if len(j.IndistinguishableProtein) > 0 {
				for _, k := range j.IndistinguishableProtein {
					ptid.IndistinguishableProtein = append(ptid.IndistinguishableProtein, k.ProteinName)
				}
			}

			for _, k := range j.Peptide {

				var pepid PeptideIonIdentification

				pepid.PeptideSequence = string(k.PeptideSequence)
				pepid.ModifiedPeptide = string(k.ModificationInfo.ModifiedPeptide)
				pepid.Charge = k.Charge
				pepid.InitialProbability = k.InitialProbability
				pepid.Weight = k.Weight
				pepid.GroupWeight = k.GroupWeight
				pepid.CalcNeutralPepMass = k.CalcNeutralPepMass
				pepid.SharedParentProteins = len(k.PeptideParentProtein)
				pepid.Modifications.Index = make(map[string]mod.Modification)
				pepid.NumberOfInstances = k.NIstances
				pepid.NumberOfEnzymaticTermini = k.NEnzymaticTermini
				pepid.Razor = -1

				if strings.EqualFold(string(k.IsNondegenerateEvidence), "Y") || strings.EqualFold(string(k.IsNondegenerateEvidence), "y") {
					pepid.IsNondegenerateEvidence = true
					pepid.IsUnique = true
				} else {
					pepid.IsNondegenerateEvidence = false
					pepid.IsUnique = false
				}

				// collect other proteins where this paptide maps to (this is different from the indistinguishable proteins list)
				for _, l := range k.PeptideParentProtein {
					pepid.PeptideParentProtein = append(pepid.PeptideParentProtein, string(l.ProteinName))
				}

				ptid.PeptideIons = append(ptid.PeptideIons, pepid)

				// get hte highest initial probability from all peptides
				if pepid.InitialProbability > ptid.TopPepProb {
					ptid.TopPepProb = pepid.InitialProbability
				}

				pepid = PeptideIonIdentification{}
			}

			peps := strings.Split(string(j.UniqueStrippedPeptides), "+")
			ptid.UniqueStrippedPeptides = peps
			proteinList = append(proteinList, ptid)
		}

		gi.Proteins = proteinList
		groups = append(groups, gi)
	}

	p.RunOptions = string(xml.ProteinSummary.ProteinSummaryHeader.ProgramDetails.ProteinProphetDetails.RunOptions)
	p.Groups = groups

	if len(groups) == 0 {
		msg.NoProteinFound(errors.New(""), "fatal")
	}

}

// PromoteProteinIDs promotes protein identifications where the reference protein
// is indistinguishable to other target proteins.
func (p *ProtXML) PromoteProteinIDs() {

	for i := range p.Groups {
		for j := range p.Groups[i].Proteins {

			var list []string
			var ref string

			if strings.HasPrefix(string(p.Groups[i].Proteins[j].ProteinName), p.DecoyTag) {
				for k := range p.Groups[i].Proteins[j].IndistinguishableProtein {
					if !strings.HasPrefix(string(p.Groups[i].Proteins[j].IndistinguishableProtein[k]), p.DecoyTag) {
						list = append(list, string(p.Groups[i].Proteins[j].IndistinguishableProtein[k]))
					}
				}
			}

			if len(list) > 0 {
				for i := range list {
					if strings.HasPrefix(list[i], "sp|") {
						ref = list[i]
						break
					} else {
						ref = list[i]
					}
				}
				p.Groups[i].Proteins[j].ProteinName = ref
			}

		}
	}

}

// MarkUniquePeptides classifies peptides as unique based on a defined threshold
func (p *ProtXML) MarkUniquePeptides(w float64) {

	for i := range p.Groups {
		for j := range p.Groups[i].Proteins {
			for k := range p.Groups[i].Proteins[j].PeptideIons {
				if p.Groups[i].Proteins[j].PeptideIons[k].Weight >= w {
					p.Groups[i].Proteins[j].PeptideIons[k].IsUnique = true
				}
			}
		}
	}

}

// Serialize converts the whle structure to a gob file
func (p *ProtIDList) Serialize() {

	b, e := msgpack.Marshal(&p)
	if e != nil {
		msg.MarshalFile(e, "fatal")
	}

	e = ioutil.WriteFile(sys.ProBin(), b, sys.FilePermission())
	if e != nil {
		msg.WriteFile(e, "fatal")
	}
}

// SerializeToTemp converts the whle structure to a gob file and puts in a specific data set folder
func (p *ProtIDList) SerializeToTemp() string {

	// // reload the meta data
	var m met.Data

	// get current directory
	dir, e := os.Getwd()
	if e != nil {
		logrus.Info("check folder permissions")
	}

	m = met.New(dir)

	eDir := os.MkdirAll(m.Temp, 0755)
	if eDir != nil {
		log.Fatal(e)
	}

	b, e := msgpack.Marshal(&p)
	if e != nil {
		msg.MarshalFile(e, "fatal")
	}

	dest := fmt.Sprintf("%s%spro.bin", m.Temp, string(filepath.Separator))

	e = ioutil.WriteFile(dest, b, sys.FilePermission())
	if e != nil {
		msg.WriteFile(e, "fatal")
	}

	return dest
}

// Restore reads philosopher results files and restore the data sctructure
func (p *ProtIDList) Restore() {

	b, e := ioutil.ReadFile(sys.ProBin())
	if e != nil {
		msg.ReadFile(e, "fatal")
	}

	e = msgpack.Unmarshal(b, &p)
	if e != nil {
		msg.DecodeMsgPck(e, "fatal")
	}

}
