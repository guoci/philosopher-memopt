package rep

import (
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"

	"philosopher/lib/cla"
	"philosopher/lib/id"
	"philosopher/lib/mod"
	"philosopher/lib/msg"
)

// AssemblePeptideReport reports consist on ion reporting
func (evi *Evidence) AssemblePeptideReport(pep id.PepIDList, decoyTag string) {

	var list PeptideEvidenceList
	var pepSeqMap = make(map[string]bool) //is this a decoy
	var pepCSMap = make(map[string][]uint8)
	var pepInt = make(map[string]float64)
	var pepProt = make(map[string]string)
	var spectra = make(map[string][]string)
	var mappedGenes = make(map[string][]string)
	var mappedProts = make(map[string][]string)
	var bestProb = make(map[string]float64)
	var pepMods = make(map[string][]mod.Modification)

	for _, i := range pep {
		if !cla.IsDecoyPSM(i, decoyTag) {
			pepSeqMap[i.Peptide] = false
		} else {
			pepSeqMap[i.Peptide] = true
		}
	}

	for _, i := range evi.PSM {

		_, ok := pepSeqMap[i.Peptide]
		if ok {

			pepCSMap[i.Peptide] = append(pepCSMap[i.Peptide], i.AssumedCharge)
			spectra[i.Peptide] = append(spectra[i.Peptide], i.Spectrum)
			pepProt[i.Peptide] = i.Protein

			if i.Intensity > pepInt[i.Peptide] {
				pepInt[i.Peptide] = i.Intensity
			}

			for j := range i.MappedProteins {
				mappedProts[i.Peptide] = append(mappedProts[i.Peptide], j)
			}

			for j := range i.MappedGenes {
				mappedGenes[i.Peptide] = append(mappedGenes[i.Peptide], j)
			}

			for _, j := range i.Modifications.Index {
				pepMods[i.Peptide] = append(pepMods[i.Peptide], j)
			}

		}

		if i.Probability > bestProb[i.Peptide] {
			bestProb[i.Peptide] = i.Probability
		}

	}

	for k, v := range pepSeqMap {

		var pep PeptideEvidence
		pep.Spectra = make(map[string]uint8)
		pep.ChargeState = make(map[uint8]uint8)
		pep.MappedGenes = make(map[string]int)
		pep.MappedProteins = make(map[string]int)
		pep.Modifications.Index = make(map[string]mod.Modification)

		pep.Sequence = k

		pep.Probability = bestProb[k]

		for _, i := range spectra[k] {
			pep.Spectra[i] = 0
		}

		for _, i := range pepCSMap[k] {
			pep.ChargeState[i] = 0
		}

		for _, i := range mappedGenes[k] {
			pep.MappedGenes[i] = 0
		}

		for _, i := range mappedProts[k] {
			pep.MappedProteins[i] = 0
		}

		d, ok := pepProt[k]
		if ok {
			pep.Protein = d
		}

		mods, ok := pepMods[pep.Sequence]
		if ok {
			for _, j := range mods {
				_, okMod := pep.Modifications.Index[j.Index]
				if !okMod {
					pep.Modifications.Index[j.Index] = j
				}
			}
		}

		// is this a decoy ?
		pep.IsDecoy = v

		list = append(list, pep)
	}

	sort.Sort(list)
	evi.Peptides = list

}

// MetaPeptideReport report consist on ion reporting
func (evi Evidence) MetaPeptideReport(workspace, brand string, channels int, hasDecoys, hasLabels bool) {

	var header string
	output := fmt.Sprintf("%s%speptide.tsv", workspace, string(filepath.Separator))

	file, e := os.Create(output)
	if e != nil {
		msg.WriteFile(errors.New("peptide output file"), "fatal")
	}
	defer file.Close()

	// building the printing set tat may or not contain decoys
	var printSet PeptideEvidenceList
	for _, i := range evi.Peptides {
		if !hasDecoys {
			if !i.IsDecoy {
				printSet = append(printSet, i)
			}
		} else {
			printSet = append(printSet, i)
		}
	}

	header = "Peptide\tPrev AA\tNext AA\tPeptide Length\tCharges\tProbability\tSpectral Count\tIntensity\tAssigned Modifications\tObserved Modifications\tProtein\tProtein ID\tEntry Name\tGene\tProtein Description\tMapped Genes\tMapped Proteins"

	if brand == "tmt" {
		switch channels {
		case 6:
			header += "\tChannel 126\tChannel 127N\tChannel 128C\tChannel 129N\tChannel 130C\tChannel 131"
		case 10:
			header += "\tChannel 126\tChannel 127N\tChannel 127C\tChannel 128N\tChannel 128C\tChannel 129N\tChannel 129C\tChannel 130N\tChannel 130C\tChannel 131N"
		case 11:
			header += "\tChannel 126\tChannel 127N\tChannel 127C\tChannel 128N\tChannel 128C\tChannel 129N\tChannel 129C\tChannel 130N\tChannel 130C\tChannel 131N\tChannel 131C"
		case 16:
			header += "\tChannel 126\tChannel 127N\tChannel 127C\tChannel 128N\tChannel 128C\tChannel 129N\tChannel 129C\tChannel 130N\tChannel 130C\tChannel 131N\tChannel 131C\tChannel 132N\tChannel 132C\tChannel 133N\tChannel 133C\tChannel 134N"
		case 18:
			header += "\tChannel 126\tChannel 127N\tChannel 127C\tChannel 128N\tChannel 128C\tChannel 129N\tChannel 129C\tChannel 130N\tChannel 130C\tChannel 131N\tChannel 131C\tChannel 132N\tChannel 132C\tChannel 133N\tChannel 133C\tChannel 134N\tChannel 134C\tChannel 135N"
		default:
			header += ""
		}
	} else if brand == "itraq" {
		switch channels {
		case 4:
			header += "\tChannel 114\tChannel 115\tChannel 116\tChannel 117"
		case 8:
			header += "\tChannel 113\tChannel 114\tChannel 115\tChannel 116\tChannel 117\tChannel 118\tChannel 119\tChannel 121"
		default:
			header += ""
		}
	}

	header += "\n"

	// verify if the structure has labels, if so, replace the original channel names by them.
	if hasLabels {

		var labels [18]string

		if len(printSet[10].Labels.Channel1.CustomName) >= 1 {
			labels[0] = printSet[10].Labels.Channel1.CustomName
			labels[1] = printSet[10].Labels.Channel2.CustomName
			labels[2] = printSet[10].Labels.Channel3.CustomName
			labels[3] = printSet[10].Labels.Channel4.CustomName
			labels[4] = printSet[10].Labels.Channel5.CustomName
			labels[5] = printSet[10].Labels.Channel6.CustomName
			labels[6] = printSet[10].Labels.Channel7.CustomName
			labels[7] = printSet[10].Labels.Channel8.CustomName
			labels[8] = printSet[10].Labels.Channel9.CustomName
			labels[9] = printSet[10].Labels.Channel10.CustomName
			labels[10] = printSet[10].Labels.Channel11.CustomName
			labels[11] = printSet[10].Labels.Channel12.CustomName
			labels[12] = printSet[10].Labels.Channel13.CustomName
			labels[13] = printSet[10].Labels.Channel14.CustomName
			labels[14] = printSet[10].Labels.Channel15.CustomName
			labels[15] = printSet[10].Labels.Channel16.CustomName
			labels[16] = printSet[10].Labels.Channel17.CustomName
			labels[17] = printSet[10].Labels.Channel18.CustomName
		}

		if labels[0] == "NA" {
			header = strings.Replace(header, "\tChannel "+printSet[10].Labels.Channel1.Name, "", -1)
		} else {
			header = strings.Replace(header, "Channel "+printSet[10].Labels.Channel1.Name, labels[0], -1)
		}

		if labels[1] == "NA" {
			header = strings.Replace(header, "\tChannel "+printSet[10].Labels.Channel2.Name, "", -1)
		} else {
			header = strings.Replace(header, "Channel "+printSet[10].Labels.Channel2.Name, labels[1], -1)
		}

		if labels[2] == "NA" {
			header = strings.Replace(header, "\tChannel "+printSet[10].Labels.Channel3.Name, "", -1)
		} else {
			header = strings.Replace(header, "Channel "+printSet[10].Labels.Channel3.Name, labels[2], -1)
		}

		if labels[3] == "NA" {
			header = strings.Replace(header, "\tChannel "+printSet[10].Labels.Channel4.Name, "", -1)
		} else {
			header = strings.Replace(header, "Channel "+printSet[10].Labels.Channel4.Name, labels[3], -1)
		}

		if labels[4] == "NA" {
			header = strings.Replace(header, "\tChannel "+printSet[10].Labels.Channel5.Name, "", -1)
		} else {
			header = strings.Replace(header, "Channel "+printSet[10].Labels.Channel5.Name, labels[4], -1)
		}

		if labels[5] == "NA" {
			header = strings.Replace(header, "\tChannel "+printSet[10].Labels.Channel6.Name, "", -1)
		} else {
			header = strings.Replace(header, "Channel "+printSet[10].Labels.Channel6.Name, labels[5], -1)
		}

		if labels[6] == "NA" {
			header = strings.Replace(header, "\tChannel "+printSet[10].Labels.Channel7.Name, "", -1)
		} else {
			header = strings.Replace(header, "Channel "+printSet[10].Labels.Channel7.Name, labels[6], -1)
		}

		if labels[7] == "NA" {
			header = strings.Replace(header, "\tChannel "+printSet[10].Labels.Channel8.Name, "", -1)
		} else {
			header = strings.Replace(header, "Channel "+printSet[10].Labels.Channel8.Name, labels[7], -1)
		}

		if labels[8] == "NA" {
			header = strings.Replace(header, "\tChannel "+printSet[10].Labels.Channel9.Name, "", -1)
		} else {
			header = strings.Replace(header, "Channel "+printSet[10].Labels.Channel9.Name, labels[8], -1)
		}

		if labels[9] == "NA" {
			header = strings.Replace(header, "\tChannel "+printSet[10].Labels.Channel10.Name, "", -1)
		} else {
			header = strings.Replace(header, "Channel "+printSet[10].Labels.Channel10.Name, labels[9], -1)
		}

		if labels[10] == "NA" {
			header = strings.Replace(header, "\tChannel "+printSet[10].Labels.Channel11.Name, "", -1)
		} else {
			header = strings.Replace(header, "Channel "+printSet[10].Labels.Channel11.Name, labels[10], -1)
		}

		if labels[11] == "NA" {
			header = strings.Replace(header, "\tChannel "+printSet[10].Labels.Channel12.Name, "", -1)
		} else {
			header = strings.Replace(header, "Channel "+printSet[10].Labels.Channel12.Name, labels[11], -1)
		}

		if labels[12] == "NA" {
			header = strings.Replace(header, "\tChannel "+printSet[10].Labels.Channel13.Name, "", -1)
		} else {
			header = strings.Replace(header, "Channel "+printSet[10].Labels.Channel13.Name, labels[12], -1)
		}

		if labels[13] == "NA" {
			header = strings.Replace(header, "\tChannel "+printSet[10].Labels.Channel14.Name, "", -1)
		} else {
			header = strings.Replace(header, "Channel "+printSet[10].Labels.Channel14.Name, labels[13], -1)
		}

		if labels[14] == "NA" {
			header = strings.Replace(header, "\tChannel "+printSet[10].Labels.Channel15.Name, "", -1)
		} else {
			header = strings.Replace(header, "Channel "+printSet[10].Labels.Channel15.Name, labels[14], -1)
		}

		if labels[15] == "NA" {
			header = strings.Replace(header, "\tChannel "+printSet[10].Labels.Channel16.Name, "", -1)
		} else {
			header = strings.Replace(header, "Channel "+printSet[10].Labels.Channel16.Name, labels[15], -1)
		}

		if labels[16] == "NA" {
			header = strings.Replace(header, "\tChannel "+printSet[10].Labels.Channel17.Name, "", -1)
		} else {
			header = strings.Replace(header, "Channel "+printSet[10].Labels.Channel17.Name, labels[16], -1)
		}

		if labels[17] == "NA" {
			header = strings.Replace(header, "\tChannel "+printSet[10].Labels.Channel18.Name, "", -1)
		} else {
			header = strings.Replace(header, "Channel "+printSet[10].Labels.Channel18.Name, labels[17], -1)
		}
	}

	_, e = io.WriteString(file, header)
	if e != nil {
		msg.WriteToFile(errors.New("cannot print PSM to file"), "fatal")
	}

	for _, i := range printSet {

		assL, obs := getModsList(i.Modifications.Index)

		var mappedProteins []string
		for j := range i.MappedProteins {
			if j != i.Protein {
				mappedProteins = append(mappedProteins, j)
			}
		}

		var cs []string
		for j := range i.ChargeState {
			cs = append(cs, strconv.Itoa(int(j)))
		}

		var mappedGenes []string
		for j := range i.MappedGenes {
			if j != i.GeneName && len(j) > 0 {
				mappedGenes = append(mappedGenes, j)
			}
		}

		sort.Strings(mappedGenes)
		sort.Strings(mappedProteins)
		sort.Strings(assL)
		sort.Strings(obs)
		sort.Strings(cs)

		line := fmt.Sprintf("%s\t%s\t%s\t%d\t%s\t%.4f\t%d\t%f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
			i.Sequence,
			i.PrevAA,
			i.NextAA,
			len(i.Sequence),
			strings.Join(cs, ", "),
			i.Probability,
			i.Spc,
			i.Intensity,
			strings.Join(assL, ", "),
			strings.Join(obs, ", "),
			i.Protein,
			i.ProteinID,
			i.EntryName,
			i.GeneName,
			i.ProteinDescription,
			strings.Join(mappedGenes, ", "),
			strings.Join(mappedProteins, ", "),
		)

		switch channels {
		case 4:

			if i.Labels.Channel1.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel1.Intensity)
			}

			if i.Labels.Channel2.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel2.Intensity)
			}

			if i.Labels.Channel3.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel3.Intensity)
			}

			if i.Labels.Channel4.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel4.Intensity)
			}

		case 6:

			if i.Labels.Channel1.CustomName != "NA" {
				line = fmt.Sprintf("\t%.4f", i.Labels.Channel1.Intensity)
			}

			if i.Labels.Channel2.CustomName != "NA" {
				line = fmt.Sprintf("\t%.4f", i.Labels.Channel2.Intensity)
			}

			if i.Labels.Channel5.CustomName != "NA" {
				line = fmt.Sprintf("\t%.4f", i.Labels.Channel5.Intensity)
			}

			if i.Labels.Channel6.CustomName != "NA" {
				line = fmt.Sprintf("\t%.4f", i.Labels.Channel6.Intensity)
			}

			if i.Labels.Channel9.CustomName != "NA" {
				line = fmt.Sprintf("\t%.4f", i.Labels.Channel9.Intensity)
			}

			if i.Labels.Channel10.CustomName != "NA" {
				line = fmt.Sprintf("\t%.4f", i.Labels.Channel10.Intensity)
			}

		case 8:

			if i.Labels.Channel1.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel1.Intensity)
			}

			if i.Labels.Channel2.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel2.Intensity)
			}

			if i.Labels.Channel3.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel3.Intensity)
			}

			if i.Labels.Channel4.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel4.Intensity)
			}

			if i.Labels.Channel5.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel5.Intensity)
			}

			if i.Labels.Channel6.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel6.Intensity)
			}

			if i.Labels.Channel7.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel7.Intensity)
			}

			if i.Labels.Channel8.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel8.Intensity)
			}

		case 10:

			if i.Labels.Channel1.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel1.Intensity)
			}

			if i.Labels.Channel2.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel2.Intensity)
			}

			if i.Labels.Channel3.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel3.Intensity)
			}

			if i.Labels.Channel4.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel4.Intensity)
			}

			if i.Labels.Channel5.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel5.Intensity)
			}

			if i.Labels.Channel6.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel6.Intensity)
			}

			if i.Labels.Channel7.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel7.Intensity)
			}

			if i.Labels.Channel8.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel8.Intensity)
			}

			if i.Labels.Channel9.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel9.Intensity)
			}

			if i.Labels.Channel10.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel10.Intensity)
			}

		case 11:

			if i.Labels.Channel1.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel1.Intensity)
			}

			if i.Labels.Channel2.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel2.Intensity)
			}

			if i.Labels.Channel3.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel3.Intensity)
			}

			if i.Labels.Channel4.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel4.Intensity)
			}

			if i.Labels.Channel5.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel5.Intensity)
			}

			if i.Labels.Channel6.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel6.Intensity)
			}

			if i.Labels.Channel7.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel7.Intensity)
			}

			if i.Labels.Channel8.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel8.Intensity)
			}

			if i.Labels.Channel9.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel9.Intensity)
			}

			if i.Labels.Channel10.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel10.Intensity)
			}

			if i.Labels.Channel11.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel11.Intensity)
			}

		case 16:

			if i.Labels.Channel1.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel1.Intensity)
			}

			if i.Labels.Channel2.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel2.Intensity)
			}

			if i.Labels.Channel3.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel3.Intensity)
			}

			if i.Labels.Channel4.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel4.Intensity)
			}

			if i.Labels.Channel5.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel5.Intensity)
			}

			if i.Labels.Channel6.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel6.Intensity)
			}

			if i.Labels.Channel7.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel7.Intensity)
			}

			if i.Labels.Channel8.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel8.Intensity)
			}

			if i.Labels.Channel9.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel9.Intensity)
			}

			if i.Labels.Channel10.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel10.Intensity)
			}

			if i.Labels.Channel11.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel11.Intensity)
			}

			if i.Labels.Channel12.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel12.Intensity)
			}

			if i.Labels.Channel13.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel13.Intensity)
			}

			if i.Labels.Channel14.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel14.Intensity)
			}

			if i.Labels.Channel15.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel15.Intensity)
			}

			if i.Labels.Channel16.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel16.Intensity)
			}

		case 18:

			if i.Labels.Channel1.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel1.Intensity)
			}

			if i.Labels.Channel2.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel2.Intensity)
			}

			if i.Labels.Channel3.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel3.Intensity)
			}

			if i.Labels.Channel4.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel4.Intensity)
			}

			if i.Labels.Channel5.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel5.Intensity)
			}

			if i.Labels.Channel6.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel6.Intensity)
			}

			if i.Labels.Channel7.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel7.Intensity)
			}

			if i.Labels.Channel8.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel8.Intensity)
			}

			if i.Labels.Channel9.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel9.Intensity)
			}

			if i.Labels.Channel10.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel10.Intensity)
			}

			if i.Labels.Channel11.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel11.Intensity)
			}

			if i.Labels.Channel12.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel12.Intensity)
			}

			if i.Labels.Channel13.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel13.Intensity)
			}

			if i.Labels.Channel14.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel14.Intensity)
			}

			if i.Labels.Channel15.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel15.Intensity)
			}

			if i.Labels.Channel16.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel16.Intensity)
			}

			if i.Labels.Channel17.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel17.Intensity)
			}

			if i.Labels.Channel18.CustomName != "NA" {
				line = fmt.Sprintf("%s\t%.4f", line, i.Labels.Channel18.Intensity)
			}

		default:
			header += ""
		}

		line += "\n"

		_, e = io.WriteString(file, line)
		if e != nil {
			msg.WriteToFile(errors.New("cannot print Peptides to file"), "fatal")
		}
	}
}
