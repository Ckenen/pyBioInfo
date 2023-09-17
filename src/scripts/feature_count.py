#!/usr/scripts/env python

import sys
from optparse import OptionParser
from pyBioInfo.IO import BamFile, BedFile
from pyBioInfo.Utils import ShiftLoader

MODE_COINCIDE = 1
MODE_OVERLAP = 2


def run(options, args):
    infile1, infile2, outfile = args
    with BamFile(infile1) as handle1, BedFile(infile2) as handle2, open(outfile, "w+") as fw:
        fw.write("TranscriptID\tLength\tCount\n")

        mode = MODE_COINCIDE
        if options.mode == "coincide":
            mode = MODE_COINCIDE
        elif options.mode == "overlap":
            mode = MODE_OVERLAP
        else:
            raise ValueError()

        strand = None
        if options.strand is None:
            strand = None
        if options.strand == "forward":
            strand = True
        elif options.strand == "reverse":
            strand = False
        else:
            raise ValueError()

        paired = False
        loader = None
        if options.type == "paired":
            paired = True
            loader = ShiftLoader(handle1.fragments())
        elif options.type == "single":
            paired = False
            loader = ShiftLoader(handle1.alignments())
        else:
            raise ValueError()

        for transcript in handle2:
            count = 0
            length = len(transcript)
            for item in loader.fetch(transcript.chrom, transcript.start, transcript.end):
                if strand is not None:
                    if strand:
                        if item.strand != transcript.strand:
                            continue
                    else:
                        if item.strand == transcript.strand:
                            continue
                if paired:
                    if mode == MODE_COINCIDE:
                        if not transcript.coincide(item.mate1):
                            continue
                        if not transcript.coincide(item.mate2):
                            continue
                    elif mode == MODE_OVERLAP:
                        if not transcript.overlap(item):
                            continue
                    else:
                        raise ValueError()
                else:
                    if mode == 1:
                        if not transcript.coincide(item):
                            continue
                    elif mode == 2:
                        if not transcript.overlap(item):
                            continue
                    else:
                        raise ValueError()
                count += 1
            fw.write("\t".join(map(str, [transcript.name, length, count])) + "\n")


def parameters():
    parser = OptionParser(usage="%prog [options] infile.fam input.bed outfile.tsv")
    parser.add_option("-s", "--strand", dest="strand", default=None)
    parser.add_option("-m", "--mode", dest="mode", default="coincide")
    parser.add_option("-t", "--type", dest="type", default="paired")
    options, args = parser.parse_args()
    if len(args) != 3:
        parser.print_help()
        exit(2)
    return options, args


def main():
    options, args = parameters()
    run(options, args)


if __name__ == '__main__':
    main()
