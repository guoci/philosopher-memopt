package rep

import (
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"philosopher/lib/msg"

	"philosopher/lib/bio"
	"philosopher/lib/cla"
	"philosopher/lib/id"
	"philosopher/lib/mod"
	"philosopher/lib/uti"
)

// AssembleIonReport reports consist on ion reporting
func (evi *Evidence) AssembleIonReport(ion id.PepIDList, decoyTag string) {

	var list IonEvidenceList
	var psmPtMap = make(map[string][]string)
	var psmIonMap = make(map[string][]string)
	var bestProb = make(map[string]float64)

	var ionMods = make(map[string][]mod.Modification)

	// collapse all psm to protein based on Peptide-level identifications
	for _, i := range evi.PSM {

		psmIonMap[i.IonForm] = append(psmIonMap[i.IonForm], i.Spectrum)
		psmPtMap[i.Spectrum] = append(psmPtMap[i.Spectrum], i.Protein)

		if i.Probability > bestProb[i.IonForm] {
			bestProb[i.IonForm] = i.Probability
		}

		for j := range i.MappedProteins {
			psmPtMap[i.IonForm] = append(psmPtMap[i.IonForm], j)
		}

		for _, j := range i.Modifications.Index {
			ionMods[i.IonForm] = append(ionMods[i.IonForm], j)
		}

	}

	for _, i := range ion {
		var pr IonEvidence

		pr.IonForm = fmt.Sprintf("%s#%d#%.4f", i.Peptide, i.AssumedCharge, i.CalcNeutralPepMass)

		pr.Spectra = make(map[string]int)
		pr.MappedGenes = make(map[string]int)
		pr.MappedProteins = make(map[string]int)
		pr.Modifications.Index = make(map[string]mod.Modification)

		v, ok := psmIonMap[pr.IonForm]
		if ok {
			for _, j := range v {
				pr.Spectra[j]++
			}
		}

		pr.Sequence = i.Peptide
		pr.ModifiedSequence = i.ModifiedPeptide
		pr.MZ = uti.Round(((i.CalcNeutralPepMass + (float64(i.AssumedCharge) * bio.Proton)) / float64(i.AssumedCharge)), 5, 4)
		pr.ChargeState = i.AssumedCharge
		pr.PeptideMass = i.CalcNeutralPepMass
		pr.PrecursorNeutralMass = i.PrecursorNeutralMass
		pr.Expectation = i.Expectation
		pr.NumberOfEnzymaticTermini = i.NumberOfEnzymaticTermini
		pr.Protein = i.Protein
		pr.MappedProteins[i.Protein] = 0
		pr.Modifications = i.Modifications
		pr.Probability = bestProb[pr.IonForm]

		// get the mapped proteins
		for _, j := range psmPtMap[pr.IonForm] {
			pr.MappedProteins[j] = 0
		}

		mods, ok := ionMods[pr.IonForm]
		if ok {
			for _, j := range mods {
				_, okMod := pr.Modifications.Index[j.Index]
				if !okMod {
					pr.Modifications.Index[j.Index] = j
				}
			}
		}

		// is this bservation a decoy ?
		if cla.IsDecoyPSM(i, decoyTag) {
			pr.IsDecoy = true
		}

		list = append(list, pr)
	}

	sort.Sort(list)
	evi.Ions = list

}

// MetaIonReport reports consist on ion reporting
func (evi Evidence) MetaIonReport(workspace, brand string, channels int, hasDecoys, hasLabels bool) {

	var header string
	output := fmt.Sprintf("%s%sion.tsv", workspace, string(filepath.Separator))

	file, e := os.Create(output)
	if e != nil {
		msg.WriteFile(errors.New("peptide ion output file"), "fatal")
	}
	defer file.Close()

	// building the printing set tat may or not contain decoys
	var printSet IonEvidenceList
	for _, i := range evi.Ions {
		// This inclusion is necessary to avoid unexistent observations from being included after using the filter --mods options
		if i.Probability > 0 {
			if !hasDecoys {
				if !i.IsDecoy {
					printSet = append(printSet, i)
				}
			} else {
				printSet = append(printSet, i)
			}
		}
	}

	header = "Peptide Sequence\tModified Sequence\tPrev AA\tNext AA\tPeptide Length\tM/Z\tCharge\tObserved Mass\tProbability\tExpectation\tSpectral Count\tIntensity\tAssigned Modifications\tObserved Modifications\tProtein\tProtein ID\tEntry Name\tGene\tProtein Description\tMapped Genes\tMapped Proteins"

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
		msg.WriteToFile(errors.New("cannot print Ion to file"), "fatal")
	}

	for _, i := range printSet {

		assL, obs := getModsList(i.Modifications.Index)

		var mappedProteins []string
		for j := range i.MappedProteins {
			if j != i.Protein {
				mappedProteins = append(mappedProteins, j)
			}
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

		line := fmt.Sprintf("%s\t%s\t%s\t%s\t%d\t%.4f\t%d\t%.4f\t%.4f\t%.14f\t%d\t%.4f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
			i.Sequence,
			i.ModifiedSequence,
			i.PrevAA,
			i.NextAA,
			len(i.Sequence),
			i.MZ,
			i.ChargeState,
			i.PeptideMass,
			i.Probability,
			i.Expectation,
			len(i.Spectra),
			i.Intensity,
			strings.Join(assL, ", "),
			strings.Join(obs, ", "),
			i.Protein,
			i.ProteinID,
			i.EntryName,
			i.GeneName,
			i.ProteinDescription,
			strings.Join(mappedGenes, ","),
			strings.Join(mappedProteins, ","),
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
			msg.WriteToFile(errors.New("cannot print Ions to file"), "fatal")
		}
	}

}
