#!/usr/bin/env python
import os
import argparse

class FastqTools(object):
    @staticmethod
    def dedup(args):
        threads = args.threads
        paths = args.in_and_out
        infile1 = infile2 = outfile1 = outfile2 = None
        if len(paths) == 2:
            infile1, outfile1 = paths
        elif len(paths) == 4:
            infile1, infile2, outfile1, outfile2 = paths
        else:
            exit(1)
        
        tmpdir = outfile1 + ".tmp"
        os.makedirs(tmpdir, exist_ok=True)
        
        # fastq -> txt
        
        # sort txt
        
        # txt -> fastq
        
        # clean
        
        exit(1)
    
        
    @staticmethod
    def fastq2row(args):
        pass
    
    
    @staticmethod
    def row2fastq(args):
        pass
    

def main():
    parser = argparse.ArgumentParser(description="A toolkit to process FAM file")

    subparsers = parser.add_subparsers(title="Available subcommands",dest="command")

    dedup_parser = subparsers.add_parser("dedup", help="Dedup FASTQ files", description="Dedup FASTQ files.")
    dedup_parser.add_argument("in_and_out", nargs="+", help="Input and output FASTQ files")
    dedup_parser.add_argument("-@", "--threads", type=int, default=4, help="Number of threads (default: 4)")
    dedup_parser.set_defaults(func=FastqTools.dedup)

    # split_parser = subparsers.add_parser("split", help="Split FAM file by seqname", description="Split FAM file by seqname.")
    # split_parser.add_argument("in.fam", help="Input FAM file")
    # split_parser.add_argument("outdir", help="Output directory")
    # split_parser.set_defaults(func=FamTools.split)

    # merge_parser = subparsers.add_parser("merge", help="Merge splitted FAM files", description="Merge splitted FAM files.")
    # merge_parser.add_argument("input", help="Input SAM/BAM file")
    # merge_parser.set_defaults(func=FamTools.merge)

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return

    args.func(args)    

if __name__ == "__main__":
    main()