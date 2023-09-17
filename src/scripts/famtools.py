#!/usr/bin/env python
import sys
import os
import optparse
import pysam
from pyBioInfo.IO.File import BedFile, FamFile, VcfFile
from pyBioInfo.Utils import ShiftLoader, BundleBuilder

class FamTools(object):
    @classmethod
    def build(cls, args):
        usage = "%prog input.bam output.fam"
        parser = optparse.OptionParser(usage=usage)
        options, args = parser.parse_args(args)
        infile, outfile = args
        
        with FamFile(infile) as f, FamFile(outfile, "wb", f) as fw:
            for obj in f:
                fw.write(obj)

    @classmethod
    def split(cls, args):
        usage = "%prog [options] input.fam outdir"
        parser = optparse.OptionParser(usage=usage)
        parser.add_option("-a", "--all", dest="all", action="store_true", default=False, help="output all references.")
        options, args = parser.parse_args(args)
        
        infile, outdir = args
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        if options.all:
            fws = dict()
            with FamFile(infile) as f:
                for chrom in f.handle.references:
                    fws[chrom] = FamFile(outdir + "/%s.fam" % chrom, "wb", f)
                for frag in f:
                    fws[frag.chrom].write(frag)
                for fw in fws.values():
                    fw.close()
        else:
            chrom = None
            fw = None
            with FamFile(infile) as f:
                for frag in f:
                    if chrom is None or frag.chrom != chrom:
                        if chrom is not None:
                            fw.close()
                        chrom = frag.chrom
                        fw = FamFile(outdir + "/%s.fam" % chrom, "wb", f)
                    fw.write(frag)
            if fw is not None:
                fw.close()
            
    @classmethod
    def merge(cls, args):
        usage = "%prog input1.fam input2.fam ... output.fam"
        parser = optparse.OptionParser(usage=usage)
        options, args = parser.parse_args(args)
        assert len(args) >= 2
        infiles = args[:-1]
        outfile = args[-1]
        fw = None
        for infile in infiles:
            with pysam.AlignmentFile(infile) as f:
                if fw is None:
                    fw = pysam.AlignmentFile(outfile, "wb", f)
                for segment in f:
                    fw.write(segment)
        fw.close()
                
        
    @classmethod
    def maskevent(cls, args):
        usage = "%prog input.fam snps.bed/vcf output.fam"
        parser = optparse.OptionParser(usage=usage)
        options, args = parser.parse_args(args)
        infile1, infile2, outfile = args
        
        with FamFile(infile1) as f, FamFile(outfile, "wb", f) as fw:
            if infile2.endswith(".bed") or infile2.endswith(".bed.gz"):
                snps = BedFile(infile2)
            elif infile2.endswith(".vcf") or infile2.endswith(".vcf.gz"):
                snps = VcfFile(infile2)
            loader = ShiftLoader(snps)
            for frag in f:
                positions = set([snp.start for snp in loader.fetch(obj=frag)])
                for align in [frag.mate1, frag.mate2]:
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
                fw.write(frag)
            
            # for bundle in BundleBuilder(f, keep=True):
            #     chrom = bundle.chrom
            #     start = bundle.start_min
            #     end = bundle.end_max
            #     positions = set([snp.start for snp in loader.fetch(chrom=chrom, start=start, end=end)])
            #     for frag in bundle.data:
            #         for align in [frag.mate1, frag.mate2]:
            #             segment = align.segment
            #             events = []
            #             for item in segment.get_tag("ME").split(";"):
            #                 if item == "":
            #                     continue
            #                 e = item.split(",")
            #                 if e[1] == "-" or e[2] == "-":
            #                     continue
            #                 if int(e[0]) in positions:
            #                     continue
            #                 events.append(item)
            #             ce = ";".join(events)
            #             segment.set_tag("CE", ce)
            #         fw.write(frag)
            

def main():
    usage = """Usage:
    famtools.py command [args]

Commands:
    -h, --help  show this help.
    build       build fragmented file (FAM format).
    split       split fam by reference.
    merge       merge fam.
    maskevent   mask SNPs mismatch events. (CE tag)
"""

    if len(sys.argv) < 2:
        sys.stdout.write(usage + "\n")
        exit(1)

    command = sys.argv[1]
    args = sys.argv[2:] 
    if command == "-h":
        sys.stdout.write(usage + "\n")
        exit(1)
    elif command == "build":
        FamTools.build(args)
    elif command == "split":
        FamTools.split(args)
    elif command == "merge":
        FamTools.merge(args)
    elif command == "maskevent":
        FamTools.maskevent(args)
    else:
        sys.stderr.write("Unknown command %s\n" % command)
        exit(1)


if __name__ == "__main__":
    main()
