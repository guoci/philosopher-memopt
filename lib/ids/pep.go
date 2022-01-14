package ids

import (
	"bufio"
	"log"
	"os"
	"philosopher/lib/ptm"
	"regexp"
)

// PepXML high-level data
type PepXML struct {
	FileName     string
	SearchEngine string
	Summary      string
	PSMs         []PeptideSpectrumMatch
}

// PeptideIdentification represents a single spectrum query from a pep.xml
type PeptideSpectrumMatch struct {
	AssumedCharge                        string
	HitRank                              string
	MissedCleavages                      string
	NumberTolTerm                        string
	NumberTotalProteins                  string
	Scan                                 string
	NumberofMissedCleavages              string
	Spectrum                             string
	SpectrumFile                         string
	Peptide                              string
	ModifiedPeptide                      string
	Protein                              string
	ProteinDescription                   string
	CompesationVoltage                   string
	PrevAA                               string
	NextAA                               string
	LocalizationRange                    string
	MSFragerLocalization                 string
	MSFraggerLocalizationScoreWithPTM    string
	MSFraggerLocalizationScoreWithoutPTM string
	UncalibratedPrecursorNeutralMass     string
	PrecursorNeutralMass                 string
	PrecursorExpMass                     string
	RetentionTime                        string
	CalcNeutralPepMass                   string
	Massdiff                             string
	Probability                          string
	Expectation                          string
	DiscriminantValue                    string
	Intensity                            string
	IonMobility                          string
	Modifications                        ptm.Modifications
	//LocalizedPTMMassDiff                 map[string]string
	//AlternativeProteins                  map[string]int
	//LocalizedPTMSites                    map[string]int
}

// PepIDList is a list of PeptideSpectrumMatch
type PepIDList []PeptideSpectrumMatch

// Len function for Sort
func (p PepIDList) Len() int {
	return len(p)
}

// Less function for Sort
func (p PepIDList) Less(i, j int) bool {
	return p[i].Probability > p[j].Probability
}

// Swap function for Sort
func (p PepIDList) Swap(i, j int) {
	p[i], p[j] = p[j], p[i]
}

func (p *PepXML) Read(f string) {

	file, err := os.Open(f)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	scanner.Split(bufio.ScanLines)

	regProphet := regexp.MustCompile(`^\<peptideprophet_summary`)

	// spectrum query attributes
	regStartSQ := regexp.MustCompile(`^\<spectrum_query`)
	regEndSQ := regexp.MustCompile(`^\<\/spectrum_query>`)
	regScan := regexp.MustCompile(`start_scan\=\"(.+?)\".+`)
	regCharge := regexp.MustCompile(`assumed_charge\=\"(\d)\".+`)
	regSpectrum := regexp.MustCompile(`spectrum\=\"(.+?)\".+`)
	regPrecursorNM := regexp.MustCompile(`precursor_neutral_mass\=\"(.+?)\".+`)
	regRT := regexp.MustCompile(`retention_time_sec\=\"(.+?)\".+`)

	// search hit atributes
	regStartSH := regexp.MustCompile(`^\<search_hit`)
	regEndSH := regexp.MustCompile(`^\<\/search_hit>`)
	regPeptide := regexp.MustCompile(`peptide\=\"(.+?)\".+`)
	regMassDiff := regexp.MustCompile(`massdiff\=\"(.+?)\".+`)
	regCalcNPM := regexp.MustCompile(`calc_neutral_pep_mass\=\"(.+?)\".+`)
	regPepNextAA := regexp.MustCompile(`peptide_next_aa\=\"(\w)\".+`)
	regPepPrevAA := regexp.MustCompile(`peptide_prev_aa\=\"(\w)\".+`)
	regNumMC := regexp.MustCompile(`num_missed_cleavages\=\"(\d+)\".+`)
	regNumTT := regexp.MustCompile(`num_tol_term\=\"(.+?)\".+`)
	regNumTP := regexp.MustCompile(`num_tot_proteins\=\"(\d+)\".+`)
	regHitRank := regexp.MustCompile(`hit_rank\=\"(\d)\".+`)
	regProtein := regexp.MustCompile(`protein\=\"(.+?)\".+`)
	regProteinDescr := regexp.MustCompile(`protein_descr\=\"(.+?)\".+`)

	// Scores
	regExpect := regexp.MustCompile(`^\<search_score name\=\"expect\" value=\"(.+?)\"\/\>`)
	regProbability := regexp.MustCompile(`^\<peptideprophet_result probability\=\"(.+?)\"`)

	// Modifications
	//regAAMod := regexp.MustCompile(`\<aminoacid_modification aminoacid\=\"(\w)\" massdiff\=\"(.+?)\" mass\=\"(.+?)\" variable\=\"(\w)\"/>`)
	regStartPTM := regexp.MustCompile(`\<modification_info modified_peptide\=\"(.+?)\"\>`)
	regEndPTM := regexp.MustCompile(`\<\/modification_info\>`)
	//regModAA := regexp.MustCompile(`\<mod_aminoacid_mass mass\=\"(.+?)\" position=\"(\d+)\"/>`)

	var flag = 0
	var psm PeptideSpectrumMatch

	for scanner.Scan() {

		if regProphet.MatchString(scanner.Text()) {
			p.Summary = "PeptideProphet"
		}

		if regEndSQ.MatchString((scanner.Text())) {
			p.PSMs = append(p.PSMs, psm)
			psm = PeptideSpectrumMatch{}
		}

		if regStartSQ.MatchString((scanner.Text())) {

			flag = 1

			psm = PeptideSpectrumMatch{}

			// Scan number
			psm.Scan = regScan.FindStringSubmatch(scanner.Text())[1]

			// Assumed charge
			psm.AssumedCharge = regCharge.FindStringSubmatch(scanner.Text())[1]

			// Spectrum name
			psm.Spectrum = regSpectrum.FindStringSubmatch(scanner.Text())[1]

			// Precursor Neutral Mass
			psm.PrecursorNeutralMass = regPrecursorNM.FindStringSubmatch(scanner.Text())[1]

			// Precursor Neutral Mass
			psm.RetentionTime = regRT.FindStringSubmatch(scanner.Text())[1]

		}

		if flag == 1 && regEndSH.MatchString((scanner.Text())) {
			flag = 0
		}

		if flag == 1 && regStartSH.MatchString((scanner.Text())) {

			// Peptide
			psm.Peptide = regPeptide.FindStringSubmatch(scanner.Text())[1]

			// Massdiff
			psm.Massdiff = regMassDiff.FindStringSubmatch(scanner.Text())[1]

			// Calculated Neutral Peptide Mass
			psm.CalcNeutralPepMass = regCalcNPM.FindStringSubmatch(scanner.Text())[1]

			// Peptide Next AA
			psm.NextAA = regPepNextAA.FindStringSubmatch(scanner.Text())[1]

			// Peptide Previous AA
			psm.PrevAA = regPepPrevAA.FindStringSubmatch(scanner.Text())[1]

			// Number of Missed Cleavages
			psm.NumberofMissedCleavages = regNumMC.FindStringSubmatch(scanner.Text())[1]

			// Number of Tol. Terms
			psm.NumberTolTerm = regNumTT.FindStringSubmatch(scanner.Text())[1]

			// Number of Total Proteins
			psm.NumberTotalProteins = regNumTP.FindStringSubmatch(scanner.Text())[1]

			// Hit Rank
			psm.HitRank = regHitRank.FindStringSubmatch(scanner.Text())[1]

			// Protein
			psm.Protein = regProtein.FindStringSubmatch(scanner.Text())[1]

			// Protein Description
			psm.ProteinDescription = regProteinDescr.FindStringSubmatch(scanner.Text())[1]

		}

		// Expectation value
		if flag == 1 && regExpect.MatchString((scanner.Text())) {
			psm.Expectation = regExpect.FindStringSubmatch(scanner.Text())[1]
		}

		// Peptideprophet probability
		if flag == 1 && regProbability.MatchString((scanner.Text())) {
			psm.Probability = regProbability.FindStringSubmatch(scanner.Text())[1]
		}

		// Modifications
		// if regAAMod.MatchString((scanner.Text())) {

		// }

		if flag == 2 && regEndPTM.MatchString((scanner.Text())) {
			flag = 1
		}

		if flag == 1 && regStartPTM.MatchString((scanner.Text())) {
			flag = 2

			//psm.Modifications = make(map[string]ptm.Modification)

			//psm.ModifiedPeptide = regStartPTM.FindStringSubmatch(scanner.Text())[1]
		}

	}
}
