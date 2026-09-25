#!/usr/bin/env python
import sys
import os
import optparse
import argparse
import pysam
from pyBioInfo.IO.File import BedFile, FamFile, VcfFile
from pyBioInfo.Utils import ShiftLoader, BundleBuilder

class FamTools(object):
    @staticmethod
    def build(args):
        infile = args.__dict__["in.bam"]
        outfile = args.__dict__["out.fam"]
        
        
        
        # usage = "%prog input.bam output.fam"
        # parser = optparse.OptionParser(usage=usage)
        # options, args = parser.parse_args(args)
        # infile, outfile = args
        
        # with FamFile(infile) as f, FamFile(outfile, "wb", f) as fw:
        #     for obj in f:
        #         fw.write(obj)

    @staticmethod
    def split(args):
        infam = args.__dict__["in.fam"]
        outdir = args.outdir
        
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
            
    @staticmethod
    def merge(args):
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
        
        
    @staticmethod
    def maskevent(args):
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
    parser = argparse.ArgumentParser(description="A toolkit to process FAM file")

    subparsers = parser.add_subparsers(title="Available subcommands",dest="command")

    build_parser = subparsers.add_parser("build", help="Build FAM file", description="Build FAM file.")
    build_parser.add_argument("in.bam", help="Input BAM file")
    build_parser.add_argument("out.fam", help="Output FAM file")
    build_parser.add_argument("-@", "--threads", type=int, default=4, help="Number of threads (default: 4)")
    build_parser.set_defaults(func=FamTools.build)

    split_parser = subparsers.add_parser("split", help="Split FAM file by seqname", description="Split FAM file by seqname.")
    split_parser.add_argument("in.fam", help="Input FAM file")
    split_parser.add_argument("outdir", help="Output directory")
    split_parser.set_defaults(func=FamTools.split)

    merge_parser = subparsers.add_parser("merge", help="Merge splitted FAM files", description="Merge splitted FAM files.")
    merge_parser.add_argument("input", help="Input SAM/BAM file")
    merge_parser.set_defaults(func=FamTools.merge)

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return

    args.func(args)    


if __name__ == "__main__":
    main()
