#!/usr/bin/env python
import argparse
from typing import cast
from pyBioInfo.Tools import BamTools

def main():
    parser = argparse.ArgumentParser(description="A toolkit to process FAM file")
    subparsers = parser.add_subparsers(title="Available subcommands", dest="command")

    sub_parser = subparsers.add_parser("bam-to-bed", help="Convert BAM to BED", description="Convert BAM to BED")
    _ = sub_parser.add_argument("-i", "--input", dest="input", required=True, help="Input BAM file")
    _ = sub_parser.add_argument("-o", "--output", dest="output", required=True, help="Output BED file")
    _ = sub_parser.add_argument("-m", "--merge-fragment", type=bool, default=False, help="Merge fragment")
    _ = sub_parser.add_argument("-f", "--fill-fragment-gap", type=bool, default=False, help="Fill gap between paired alignments")
    _ = sub_parser.add_argument("-b", "--build-index", type=bool, default=False, help="Build index")
    
    sub_parser = subparsers.add_parser("extract-genome-sizes", help="Extract genome sizes from BAM file", description="Extract genome sizes from BAM file")
    _ = sub_parser.add_argument("-i", "--input", dest="input", required=True, help="Input BAM file")
    _ = sub_parser.add_argument("-o", "--output", dest="output", required=True, help="Output genome sizes file")

    sub_parser = subparsers.add_parser("bam-to-bigwig", help="Convert BAM to BIGWIG", description="Convert BAM to BIGWIG")
    _ = sub_parser.add_argument("-i", "--input", help="Input BAM file")
    _ = sub_parser.add_argument("-o", "--output1", help="Output BIGWIG file")
    _ = sub_parser.add_argument("-p", "--output2", help="Output BIGWIG file (For stranded)")
    _ = sub_parser.add_argument("-m", "--merge-fragment", type=bool, default=False, help="Merge fragment")
    _ = sub_parser.add_argument("-f", "--fill-fragment-gap", type=bool, default=False, help="Fill gap between paired alignments")
    _ = sub_parser.add_argument("-s", "--scale", type=float, default=1, help="Scale factor")
    _ = sub_parser.add_argument("-S", "--stranded", type=bool, default=False, help="Seperate by strand")  


    sub_parser = subparsers.add_parser("add-nh-hi-tag", help="Add NH and HI tag to BAM file", description="Add NH and HI tag to BAM file")
    _ = sub_parser.add_argument("-i", "--input", help="Input BAM file")
    _ = sub_parser.add_argument("-o", "--output", help="Output BAM file")
    _ = sub_parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads")
    _ = sub_parser.add_argument("-p", "--proper-pair", type=bool, default=False, help="Only output proper pairs")
    _ = sub_parser.add_argument("-s", "--input-is-sorted", type=bool, default=False, help="Input BAM is sorted by read name. Otherwise it will be sorted by read name first")
    _ = sub_parser.add_argument("-S", "--sort-output", type=bool, default=False, help="Sort the output BAM by coordinate. Otherwise the output BAM is sorted by read name")

    sub_parser = subparsers.add_parser("stat-chrom-reads", help="Stat chrom reads", description="Stat chrom reads")
    _ = sub_parser.add_argument("-i", "--input", help="Input BAM file")
    _ = sub_parser.add_argument("-o", "--output", help="Output TSV/CSV file")
    _ = sub_parser.add_argument("-p", "--primary-only", type=bool, default=False, help="Only count primary alignments")
    _ = sub_parser.add_argument("-P", "--proper-pair-only", type=bool, default=False, help="Only count proper-pair alignments")
    _ = sub_parser.add_argument("-u", "--uniq-mapped-only", type=bool, default=False, help="Only count uniq-mapped alignments")
    # _ = sub_parser.add_argument("-d", "--ignore-pcr-duplicate", type=bool, default=False, help="Ignore PCR duplicate")


    args = parser.parse_args()
    command = cast(str | None, args.command)
    if command is None:
        parser.print_help()
        return
    dispatch = {
        "bam-to-bed": BamTools.bam_to_bed,
        "bam-to-bigwig": BamTools.bam_to_bigwig,
        "extract-genome-sizes": BamTools.extract_genome_sizes,
        "add-nh-hi-tag": BamTools.add_nh_hi_tag,
        "stat-chrom-reads": BamTools.stat_chrom_reads,
    }
    dispatch[command](args)


if __name__ == "__main__":
    main()
