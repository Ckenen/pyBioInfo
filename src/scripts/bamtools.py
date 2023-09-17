#!/usr/bin/env python
import sys
from optparse import OptionParser
import pysam
from pyBioInfo.IO.File import BamFile, BedFile, VcfFile
from pyBioInfo.Utils import SegmentTools, ShiftLoader


class BamTools(object):
    @classmethod
    def extract(cls, args):
        usage = "bamtools.py extract [options] input.bam output.bam"
        parser = OptionParser(usage=usage)
        parser.add_option("-c", "--chrom", dest="chrom", type="string")
        parser.add_option("-s", "--start", dest="start",
                          default=None, type="int")
        parser.add_option("-e", "--end", dest="end", default=None, type="int")
        options, args = parser.parse_args(args)

        infile, outfile = args
        chrom = options.chrom
        start = options.start
        end = options.end
        assert chrom is not None

        with pysam.AlignmentFile(infile) as f, pysam.AlignmentFile(outfile, "wb", template=f) as fw:
            for segment in f.fetch(contig=chrom, start=start, end=end):
                fw.write(segment)

    @classmethod
    def getevent(cls, args):
        usage = "bamtools.py getevent input.bam output.bam"
        parser = OptionParser(usage=usage)
        options, args = parser.parse_args(args)
        infile, outfile = args
        with pysam.AlignmentFile(infile) as f, pysam.AlignmentFile(outfile, "wb", f) as fw:
            for segment in f:
                if segment.is_unmapped or segment.is_secondary or segment.is_supplementary:
                    continue
                bases = SegmentTools.get_reference_sequence_base(segment)
                rc = ";".join(["%s,%d" % (k, bases[k]) for k in sorted(bases.keys())])
                items = []
                for e in SegmentTools.get_events(segment):
                    if isinstance(e[3], list):
                        e[3] = "/".join(map(str, e[3]))
                    items.append(",".join(map(str, e)))
                me = ";".join(items)
                segment.set_tag("RC", rc)
                segment.set_tag("ME", me)
                fw.write(segment)
                
    @classmethod
    def maskevent(cls, args):
        usage = "bamtools.py maskevent input.bam snps.bed output.bam"
        parser = OptionParser(usage=usage)
        options, args = parser.parse_args(args)
        infile1, infile2, outfile = args
        
        with BamFile(infile1) as f, BamFile(outfile, "wb", f) as fw:
            if infile2.endswith(".bed") or infile2.endswith(".bed.gz"):
                snps = BedFile(infile2)
            elif infile2.endswith(".vcf") or infile2.endswith(".vcf.gz"):
                snps = VcfFile(infile2)
            loader = ShiftLoader(snps)
            for align in f:
                positions = set([snp.start for snp in loader.fetch(obj=align)])
                segment = align.segment
                events = []
                for item in segment.get_tag("ME").split(";"):
                    if item == "":
                        continue
                    e = item.split(",")
                    if e[1] == "-" or e[2] == "-":
                        continue
                    if int(e[0]) in positions:
                        continue
                    events.append(item)
                ce = ";".join(events)
                segment.set_tag("CE", ce)
                fw.write(align)


def main():
    usage = """Usage:
    bamtools.py command [args]

Commands:
    extract       extract sub-bam.
    getevent      get mismatch events. (ME tag and RC tag)
    maskevent     mask SNPs mismatch events. (CE tag)
    """

    if len(sys.argv) < 2:
        print(usage)
        exit(1)

    command = sys.argv[1]
    args = sys.argv[2:]

    if command == "extract":
        BamTools.extract(args)
    elif command == "getevent":
        BamTools.getevent(args)
    elif command == "maskevent":
        BamTools.maskevent(args)
    else:
        print(usage)
        exit(1)


if __name__ == "__main__":
    main()
