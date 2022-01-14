package ptm

// Modifications is a collection of modifications
type Modifications map[string]Modification

// Modification is the basic attribute for each modification
type Modification struct {
	ID                string
	Name              string
	Definition        string
	Variable          string
	Position          string
	Type              string
	AminoAcid         string
	IsProteinTerminus string
	Terminus          string
	MonoIsotopicMass  float32
	AverageMass       float32
	MassDiff          float32
	IsobaricMods      map[string]float32
}
