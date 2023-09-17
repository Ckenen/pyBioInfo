#!/usr/bin/env python
import sys
import os
from pyBioInfo.IO.File import BamFile, FamFile, BedFile
from pyBioInfo.Utils import ShiftLoader


def load_snps(path):
    with BedFile(path) as f:
        for obj in f:
            yield obj


def annotate_snp_se(infile1, infile2, outfile):
    loader = ShiftLoader(load_snps(infile2))
    with BamFile(infile1) as f, BamFile(outfile, "wb", f) as fw:
        for obj in f:
            events = []
            positions = set([snp.start for snp in loader.fetch(obj=obj)])
            for item in obj.segment.get_tag("ME").split(";"):
                if item != "":
                    event = item.split(",")
                    position = int(event[0])
                    if position not in positions:
                        events.append(event)
            ce = ";".join([",".join(item) for item in events])
            obj.segment.set_tag("CE", ce)
            fw.write(obj)
            
            
def annotate_snp_pe(infile1, infile2, outfile):
    loader = ShiftLoader(load_snps(infile2))
    with FamFile(infile1) as f, FamFile(outfile, "wb", f) as fw:
        for obj in f:
            events = []
            positions = set([snp.start for snp in loader.fetch(obj=obj)])
            for item in obj.mate1.segment.get_tag("PE").split(";"):
                if item != "":
                    event = item.split(",")
                    position = int(event[0])
                    if position not in positions:
                        events.append(event)
            ce = ";".join([",".join(item) for item in events])
            obj.mate1.segment.set_tag("CE", ce)
            fw.write(obj)
            
            
def main():
    infile1 = None # identified.bam
    infile2 = None # snps.bed
    outfile = None # annotated.bam
    try:
        infile1, infile2, outfile = sys.argv[1:]
    except Exception:
        print("Usage: %s identified.bam snps.bed annotated.bam" % os.path.basename(__file__))
        exit(1)
        
    if infile1.endswith(".bam"):
        annotate_snp_se(infile1, infile2, outfile)
    elif infile1.endswith(".fam"):
        annotate_snp_pe(infile1, infile2, outfile)
    else:
        raise RuntimeError()


if __name__ == '__main__':
    main()
    
