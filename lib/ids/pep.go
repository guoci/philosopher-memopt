package ids

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"philosopher/lib/mod"
	"regexp"
	"strconv"
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
	AssumedCharge                        uint8
	HitRank                              uint8
	MissedCleavages                      uint8
	NumberTolTerm                        uint8
	NumberOfEnzymaticTermini             uint8
	NumberTotalProteins                  uint16
	Index                                uint32
	Scan                                 int
	NumberofMissedCleavages              int
	Spectrum                             string
	SpectrumFile                         string
	Peptide                              string
	ModifiedPeptide                      string
	Protein                              string
	CompesationVoltage                   string
	PrevAA                               string
	NextAA                               string
	LocalizationRange                    string
	MSFragerLocalization                 string
	MSFraggerLocalizationScoreWithPTM    string
	MSFraggerLocalizationScoreWithoutPTM string
	UncalibratedPrecursorNeutralMass     float32
	PrecursorNeutralMass                 float32
	PrecursorExpMass                     float32
	RetentionTime                        float32
	CalcNeutralPepMass                   float32
	Massdiff                             float32
	Probability                          float32
	Expectation                          float32
	DiscriminantValue                    float32
	Intensity                            float32
	IonMobility                          float32
	Modifications                        mod.Modifications
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

	regProphet := regexp.MustCompile(`^\<peptideprophet\_summary`)

	regStartSQ := regexp.MustCompile(`^\<spectrum\_query`)
	regEndSQ := regexp.MustCompile(`^\<\/spectrum\_query>`)

	regScan := regexp.MustCompile(`start\_scan\=\"(.+?)\".+`)

	//var flag = 0
	for scanner.Scan() {

		var psm PeptideSpectrumMatch

		if regProphet.MatchString(scanner.Text()) {
			p.Summary = "PeptideProphet"
		}

		if regStartSQ.MatchString((scanner.Text())) {

			var psm PeptideSpectrumMatch

			scan := regScan.FindStringSubmatch(scanner.Text())
			psm.Scan, _ = strconv.Atoi(scan[1])
			fmt.Println(psm.Scan)
			scan = nil

		}

		if regEndSQ.MatchString((scanner.Text())) {
			p.PSMs = append(p.PSMs, psm)
			psm = PeptideSpectrumMatch{}
		}
	}
}
