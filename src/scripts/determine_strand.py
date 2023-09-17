#!/usr/bin/env python
import optparse
from collections import defaultdict
import pysam
from pyBioInfo.IO.File import BamFile, FamFile, BedFile
from pyBioInfo.Utils import ShiftLoader


# def load_genes(path):
#     with BedFile(path) as f:
#         for gene in f:
#             yield gene


def load_transcripts(path):
    with BedFile(path) as f:
        for transcript in f:
            yield transcript


def process_se(infile, outfile, strand, gene):
    counter = defaultdict(int)
    if strand is None:
        with BamFile(infile) as f, BamFile(outfile, "wb", f) as fw:
            loader = ShiftLoader(load_transcripts(gene))
            for obj in f:
                st = "U"
                transcripts = [t for t in loader.fetch(obj=obj)]
                strands = set([t.strand for t in transcripts])
                if len(strands) == 0:
                    st = "U"
                elif len(strands) == 1:
                    st = list(strands)[0]
                elif len(strands) == 2:
                    temp1 = set()
                    temp2 = []
                    for t in transcripts:
                        if t.start <= obj.start and t.end >= obj.end:
                            temp1.add(t.strand)
                            temp2.append(t)
                    if len(temp1) == 0:
                        st = "A"
                    elif len(temp1) == 1:
                        st = list(temp1)[0]
                    else:
                        temp3 = set()
                        for t in temp2:
                            if t.coincide(obj):
                                temp3.add(t.strand)
                        if len(temp3) == 0:
                            st = "A"
                        elif len(temp3) == 1:
                            st = list(temp3)[0]
                        else:
                            st = "A"
                else:
                    assert False
                counter[st] += 1
                if st == "+" or st == "-":
                    obj.segment.set_tag("ST", st)
                    fw.write(obj)            
    elif strand == "fr":
        with pysam.AlignmentFile(infile) as f, pysam.AlignmentFile(outfile, "wb", f) as fw:
            for segment in f:
                if segment.is_unmapped:
                    continue
                st = "-" if segment.is_reverse else "+"
                counter[st] += 1
                segment.set_tag("ST", st)
                fw.write(segment)
    elif strand == "rf":
        with pysam.AlignmentFile(infile) as f, pysam.AlignmentFile(outfile, "wb", f) as fw:
            for segment in f:
                if segment.is_unmapped:
                    continue
                st = "+" if segment.is_reverse else "-"
                counter[st] += 1
                segment.set_tag("ST", st)
                fw.write(segment)
    else:
        raise RuntimeError()
    return counter
    

def process_pe(infile, outfile, strand, gene):
    counter = defaultdict(int)
    if strand is None:
        with FamFile(infile) as f, FamFile(outfile, "wb", f) as fw:
            loader = ShiftLoader(load_transcripts(gene))
            for obj in f:
                st = "U"
                transcripts = [t for t in loader.fetch(obj=obj)]
                strands = set([t.strand for t in transcripts])
                if len(strands) == 0:
                    st = "U"
                elif len(strands) == 1:
                    st = list(strands)[0]
                elif len(strands) == 2:
                    temp1 = set()
                    temp2 = []
                    for t in transcripts:
                        if t.start <= obj.start and t.end >= obj.end:
                            temp1.add(t.strand)
                            temp2.append(t)
                    if len(temp1) == 0:
                        st = "A"
                    elif len(temp1) == 1:
                        st = list(temp1)[0]
                    else:
                        temp3 = set()
                        for t in temp2:
                            if obj.coincide_from(t):
                                temp3.add(t.strand)
                        if len(temp3) == 0:
                            st = "A"
                        elif len(temp3) == 1:
                            st = list(temp3)[0]
                        else:
                            st = "A"
                else:
                    assert False
                counter[st] += 1
                if st == "+" or st == "-":
                    obj.mate1.segment.set_tag("ST", st)
                    obj.mate2.segment.set_tag("ST", st)
                    fw.write(obj)        
    elif strand == "fr":
        with FamFile(infile) as f, FamFile(outfile, "wb", f) as fw:
            for align in f:
                st = align.mate1.strand
                counter[st] += 1
                align.mate1.segment.set_tag("ST", st)
                align.mate2.segment.set_tag("ST", st)
                fw.write(align)
    elif strand == "rf":
        with FamFile(infile) as f, FamFile(outfile, "wb", f) as fw:
            for align in f:
                st = align.mate2.strand
                counter[st] += 1
                align.mate1.segment.set_tag("ST", st)
                align.mate2.segment.set_tag("ST", st)
                fw.write(align)
    else:
        raise RuntimeError()
    return counter
    

def main():
    parser = optparse.OptionParser()
    parser.add_option("-s", "--strand", dest="strand")
    parser.add_option("-g", "--gene", dest="gene")
    options, args = parser.parse_args()
    infile, outfile = args
    
    strand = options.strand
    gene = options.gene
    
    counter = None
    if infile.endswith(".bam"):
        counter = process_se(infile, outfile, strand, gene)
    elif infile.endswith(".fam"):
        counter = process_pe(infile, outfile, strand, gene)
    else:
        assert False
        
    name = infile.split("/")[-1][:-4]
    keys = ["+", "-", "A", "U"]
    counts = [counter[k] for k in keys]
    print("Name", "Forward", "Reverse", "Ambiguous", "Unknown", sep="\t")
    print(name, counts[0], counts[1], counts[2], counts[3], sep="\t")
        
    
if __name__ == '__main__':
    main()
    
