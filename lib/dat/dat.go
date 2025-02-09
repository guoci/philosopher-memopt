// Package dat (Database)
package dat

import (
	"errors"
	"fmt"
	"io"
	"io/ioutil"
	"net/http"
	"os"
	"path"
	"path/filepath"
	"sort"
	"strings"
	"time"

	"philosopher/lib/msg"

	"philosopher/lib/fas"
	"philosopher/lib/met"
	"philosopher/lib/sys"

	"github.com/sirupsen/logrus"
	"github.com/vmihailenco/msgpack"
)

// Base main structure
type Base struct {
	FileName        string
	UniProtDB       string
	CrapDB          string
	Prefix          string
	DownloadedFiles []string
	TaDeDB          map[string]string
	Records         []Record
}

// New constructor
func New() Base {

	var self Base

	self.TaDeDB = make(map[string]string)
	self.Records = []Record{}

	return self
}

// Run is the main entry point for the databse command
func Run(m met.Data) met.Data {

	var db = New()

	if len(m.Database.ID) == 0 && (len(m.Database.Annot) == 0 || m.Database.Annot == "--contam" || m.Database.Annot == "--prefix") && (len(m.Database.Custom) == 0 || m.Database.Custom == "--contam" || m.Database.Custom == "--prefix") {
		msg.InputNotFound(errors.New("provide a protein FASTA file or Proteome ID"), "fatal")
	}

	if len(m.Database.Annot) > 0 {

		logrus.Info("Annotating the database")

		m.DB = m.Database.Annot

		db.ProcessDB(m.Database.Annot, m.Database.Tag)

		db.Serialize()

		return m
	}

	if len(m.Database.ID) < 1 && len(m.Database.Custom) < 1 {
		msg.InputNotFound(errors.New("you need to provide a taxon ID or a custom FASTA file"), "fatal")
	}

	if !m.Database.Crap {
		msg.Custom(errors.New("contaminants are not going to be added to database"), "warning")
	}

	if len(m.Database.Custom) < 1 {

		m.DB = m.Database.Custom

		dbs := strings.Split(m.Database.ID, ",")
		for _, i := range dbs {
			logrus.Info("Fetching database ", i)

			currentTime := time.Now()
			m.Database.TimeStamp = currentTime.Format("2006.01.02 15:04:05")

			db.Fetch(i, m.Temp, m.Database.Iso, m.Database.Rev)
		}

	} else {
		dbPath, _ := filepath.Abs(m.Database.Custom)
		db.UniProtDB = dbPath
		db.DownloadedFiles = append(db.DownloadedFiles, dbPath)
	}

	logrus.Info("Generating the target-decoy database")
	db.Create(m.Temp, m.Database.Add, m.Database.Enz, m.Database.Tag, m.Database.Crap, m.Database.NoD, m.Database.CrapTag)

	logrus.Info("Creating file")
	customDB := db.Save(m.Home, m.Temp, m.Database.ID, m.Database.Tag, m.Database.Rev, m.Database.Iso, m.Database.NoD, m.Database.Crap)

	db.ProcessDB(customDB, m.Database.Tag)

	logrus.Info("Processing decoys")
	db.Create(m.Temp, m.Database.Add, m.Database.Enz, m.Database.Tag, m.Database.Crap, m.Database.NoD, m.Database.CrapTag)

	logrus.Info("Creating file")
	db.Save(m.Home, m.Temp, m.Database.ID, m.Database.Tag, m.Database.Rev, m.Database.Iso, m.Database.NoD, m.Database.Crap)

	db.Prefix = m.Database.Tag

	db.Serialize()

	return m
}

// ProcessDB determines the type of sequence and sends it to the appropriate parsing function
func (d *Base) ProcessDB(file, decoyTag string) {

	fastaMap := fas.ParseFile(file)
	d.FileName = path.Base(file)

	for k, v := range fastaMap {

		class := Classify(k, decoyTag)

		if class == "uniprot" {

			db := ProcessUniProtKB(k, v, decoyTag)
			d.Records = append(d.Records, db)

		} else if class == "ncbi" {

			db := ProcessNCBI(k, v, decoyTag)
			d.Records = append(d.Records, db)

		} else if class == "ensembl" {

			db := ProcessENSEMBL(k, v, decoyTag)
			d.Records = append(d.Records, db)

		} else if class == "generic" {

			db := ProcessGeneric(k, v, decoyTag)
			d.Records = append(d.Records, db)

		} else if class == "uniref" {

			db := ProcessUniRef(k, v, decoyTag)
			d.Records = append(d.Records, db)

		} else {
			msg.ParsingFASTA(errors.New(""), "fatal")
		}
	}

}

// Fetch downloads a database file from UniProt
func (d *Base) Fetch(id, temp string, iso, rev bool) {

	var query string

	d.UniProtDB = fmt.Sprintf("%s%s%s.fas", temp, string(filepath.Separator), id)

	if rev {
		query = fmt.Sprintf("%s%s%s", "http://www.uniprot.org/uniprot/?query=reviewed:yes+AND+proteome:", id, "&format=fasta")
	} else {
		query = fmt.Sprintf("%s%s%s", "http://www.uniprot.org/uniprot/?query=proteome:", id, "&format=fasta")
	}

	if iso {
		query = fmt.Sprintf("%s&include=yes", query)
	} else {
		query = fmt.Sprintf("%s&include=no", query)
	}

	// tries to create an output file
	output, e := os.Create(d.UniProtDB)
	if e != nil {
		msg.WriteFile(errors.New("cannot create a local database file"), "fatal")
	}
	defer output.Close()

	// Tries to query data from Uniprot
	response, e := http.Get(query)
	if e != nil {
		msg.Custom(errors.New("uniProt query failed, please check your connection"), "fatal")
	}

	if response.ContentLength != -1 {
		msg.Custom(errors.New("no sequences downloaded, check your proteome ID and parameters"), "fatal")
	}
	defer response.Body.Close()

	// Tries to download data from Uniprot
	_, e = io.Copy(output, response.Body)
	if e != nil {
		msg.Custom(errors.New("UniProt download failed, please check your connection"), "fatal")
	}

	d.DownloadedFiles = append(d.DownloadedFiles, d.UniProtDB)

}

// Create processes the given fasta file and add decoy sequences
func (d *Base) Create(temp, add, enz, tag string, crap, noD, cTag bool) {

	d.TaDeDB = make(map[string]string)

	for _, i := range d.DownloadedFiles {

		dbfile, _ := filepath.Abs(i)
		db := fas.ParseFile(dbfile)

		if len(add) > 0 {
			add := fas.ParseFile(add)

			for k, v := range add {
				db[k] = v
			}
		}

		// adding contaminants to database before reversion
		// repeated entries are removed and substituted by contaminants
		if crap {

			d.Deploy(temp)

			crapMap := fas.ParseFile(d.CrapDB)

			for k, v := range crapMap {

				if cTag {
					k = "contam_" + k
				}

				split := strings.Split(k, "|")
				for i := range db {
					if strings.Contains(i, split[1]) {
						delete(db, i)
					}
				}
				db[k] = v
			}

		}

		for h, s := range db {

			th := ">" + h
			d.TaDeDB[th] = s

			if !noD {
				dh := ">" + tag + h
				d.TaDeDB[dh] = reverseSeq(s)
			}

		}

	}

}

// Deploy crap file to session folder
func (d *Base) Deploy(temp string) {

	d.CrapDB = fmt.Sprintf("%s%scrap.fas", temp, string(filepath.Separator))

	param, e1 := Asset("crap.fas")
	if e1 != nil {
		msg.WriteFile(e1, "fatal")
	}

	e2 := ioutil.WriteFile(d.CrapDB, param, sys.FilePermission())
	if e2 != nil {
		msg.WriteFile(e2, "fatal")
	}

}

// Save fasta file to disk
func (d *Base) Save(home, temp, ids, tag string, isRev, hasIso, noD, Crap bool) string {

	var base string

	if len(ids) > 0 {
		base = strings.Replace(ids, ",", "-", -1)
	} else {
		base = filepath.Base(d.UniProtDB)
	}

	t := time.Now()
	stamp := t.Format("2006-01-02")

	baseName := fmt.Sprintf("%s%s", string(filepath.Separator), stamp)

	if !noD {
		baseName = baseName + "-decoys"
	}

	if isRev {
		baseName = baseName + "-reviewed"
	}

	if hasIso {
		baseName = baseName + "-isoforms"
	}

	if Crap {
		baseName = baseName + "-contam"
	}

	workfile := fmt.Sprintf("%s%s-%s.fas", temp, baseName, base)
	outfile := fmt.Sprintf("%s%s-%s.fas", home, baseName, base)

	// create db file
	file, e := os.Create(workfile)
	if e != nil {
		msg.ReadFile(errors.New("cannot open the database file"), "fatal")
	}
	defer file.Close()

	var headers []string
	for k := range d.TaDeDB {
		headers = append(headers, k)
	}

	sort.Strings(headers)

	for _, i := range headers {
		line := i + "\n" + d.TaDeDB[i] + "\n"
		_, e = io.WriteString(file, line)
		if e != nil {
			msg.WriteFile(e, "fatal")
		}
	}

	sys.CopyFile(workfile, outfile)

	return outfile
}

// Serialize saves to disk a msgpack version of the database data structure
func (d *Base) Serialize() {

	b, e := msgpack.Marshal(&d)
	if e != nil {
		msg.MarshalFile(e, "fatal")
	}

	e = ioutil.WriteFile(sys.DBBin(), b, sys.FilePermission())
	if e != nil {
		msg.SerializeFile(e, "fatal")
	}

}

// Restore reads philosopher results files and restore the data sctructure
func (d *Base) Restore() {

	b, e := ioutil.ReadFile(sys.DBBin())
	if e != nil {
		msg.MarshalFile(e, "warning")
	}

	e = msgpack.Unmarshal(b, &d)
	if e != nil {
		msg.SerializeFile(e, "warning")
	}

}

// RestoreWithPath reads philosopher results files and restore the data sctructure
func (d *Base) RestoreWithPath(p string) {

	path := fmt.Sprintf("%s%s%s", p, string(filepath.Separator), sys.DBBin())
	path, _ = filepath.Abs(path)

	file, _ := os.Open(path)

	dec := msgpack.NewDecoder(file)
	e := dec.Decode(&d)
	if e != nil {
		msg.DecodeMsgPck(e, "fatal")
	}

}

// reverseSeq returns its argument string reversed rune-wise left to right.
func reverseSeq(s string) string {
	r := []rune(s)
	for i, j := 0, len(r)-1; i < len(r)/2; i, j = i+1, j-1 {
		r[i], r[j] = r[j], r[i]
	}
	return string(r)
}
