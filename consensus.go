package main

import (
	"flag"
	"fmt"
	"log"
	"os"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/bed"
	"github.com/biogo/biogo/seq"
	"github.com/biogo/biogo/seq/linear"
	"github.com/biogo/biogo/seq/multi"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

// command line arguments
var (
	inbam = flag.String("bam", "", "Specifies the input bam file.")
	inbai = flag.String("bai", "", "Specifies the input bai file.")
	inbed = flag.String("bed", "", "Specifies the input bed file.")
)

const thresh = 0

func main() {
	// complain if any of the command line arguments is missing
	flag.Parse()
	if *inbed == "" {
		flag.Usage()
		os.Exit(1)
	}
	if *inbam == "" {
		flag.Usage()
		os.Exit(1)
	}
	if *inbai == "" {
		flag.Usage()
		os.Exit(1)
	}

	// open the input files (bam, bai and bed)
	inBamFile, err := os.Open(*inbam)
	if err != nil {
		log.Fatalf("Failed to open input bam file: %v", err)
	}
	defer inBamFile.Close()
	// r declares a new Reader from the input bam file
	r, err := bam.NewReader(inBamFile, 0)
	r.Omit(bam.AuxTags)
	if err != nil {
		log.Fatalf("Failed to read input bam file", err)
	}
	// refs is a map of string keys to sam refs
	refs := make(map[string]*sam.Reference)
	for _, r := range r.Header().Refs() {
		refs[r.Name()] = r
	}

	inBaiFile, err := os.Open(*inbai)
	if err != nil {
		log.Fatalf("Failed to open input bai file: %v", err)
	}
	defer inBaiFile.Close()
	// read in the bai index
	idx, err := bam.ReadIndex(inBaiFile)
	if err != nil {
		log.Fatalf("Failed to read input bai file", err)
	}

	inBedFile, err := os.Open(*inbed)
	if err != nil {
		log.Fatalf("Failed to open input bed file: %v", err)
	}
	defer inBedFile.Close()
	br, err := bed.NewReader(inBedFile, 4)
	if err != nil {
		log.Fatalf("failed to open input bed reader: %v", err)
	}

	sc := featio.NewScanner(br)
	for sc.Next() {
		f := sc.Feat()
		ref := refs[f.Location().Name()]
		chunks, err := idx.Chunks(ref, f.Start(), f.End())
		if err != nil {
			log.Fatalf("could not create chunks: %v", err)
		}
		i, err := bam.NewIterator(r, chunks)
		if err != nil {
			log.Fatalf("could not create iterator for %v: %v", f, err)
		}
		ms := &multi.Multi{ColumnConsense: seq.DefaultQConsensus}
		for i.Next() {
			r := i.Record()
			if r.MapQ >= thresh {
				// Convert r into a linear.QSeq
				qseq := make([]alphabet.QLetter, len(r.Qual))
				for i, l := range r.Seq.Expand() {
					qseq[i] = alphabet.QLetter{L: alphabet.Letter(l), Q: alphabet.Qphred(r.Qual[i])}
				}
				s := linear.NewQSeq(r.Name, qseq, alphabet.DNAgapped, alphabet.Sanger)
				s.SetOffset(r.Start())

				ms.Add(s)
			}
		}
		err = i.Close()
		if err != nil {
			log.Fatalf("failed to close iterator: %v", err)
		}

		// consensus of each column
		// upper case = high quality, lower case = low quality
		c := ms.Consensus(false)
		c.Threshold = 40
		c.QFilter = seq.CaseFilter

		fmt.Printf("%60a\n", c)
	}
	err = sc.Error()
	if err != nil {
		log.Fatalf("failed filtering bam input on bed intervals: %v", err)
	}
}
