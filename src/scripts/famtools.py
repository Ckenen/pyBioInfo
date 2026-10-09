#!/usr/bin/env python
"""
A toolkit to process FAM (fragment alignment map) files.
"""
import argparse
import os
from typing import cast
from pyBioInfo.IO.File import FamFile


class FamTools(object):
    """FAM file processing subcommands."""

    @staticmethod
    def build(args: argparse.Namespace) -> None:
        """Build a FAM file from a BAM file.

        Records are expected to be pair-adjacent (e.g. query-name
        sorted); ``FamFile`` pairs every two records into a fragment.
        """
        inbam = cast(str, args.input)
        outfam = cast(str, args.output)
        f = FamFile(inbam)
        fw = FamFile(outfam, "wb", f)
        for fragment in f.fetch():
            fw.write(fragment)
        fw.close()
        f.close()

    @staticmethod
    def split(args: argparse.Namespace) -> None:
        """Split a FAM file into one FAM file per reference sequence."""
        infam = cast(str, args.infam)
        outdir = cast(str, args.outdir)
        output_all = cast(bool, args.all)
        os.makedirs(outdir, exist_ok=True)
        
        if output_all:
            # Open one writer per reference and route each fragment to it.
            f = FamFile(infam)
            if f.references is None:
                raise RuntimeError("Failed to get reference names from the FAM file!")

            handles = {
                chrom: FamFile(os.path.join(outdir, chrom + ".fam"), "wb", f)
                for chrom in f.references
            }
            try:
                for fragment in f.fetch():
                    handles[fragment.chrom].write(fragment)
            finally:
                for fw in handles.values():
                    fw.close()
            f.close()
        else:
            # Only emit files that actually contain fragments.
            chrom_current: str | None = None
            fw: FamFile | None = None
            f = FamFile(infam)
            for fragment in f.fetch():
                if chrom_current is None or fragment.chrom != chrom_current:
                    if chrom_current is not None:
                        if fw is None:
                            raise RuntimeError("fw is None!")
                        fw.close()
                    chrom_current = fragment.chrom
                    fw = FamFile(os.path.join(outdir, chrom_current + ".fam"), "wb", f)
                if fw is None:
                    raise RuntimeError("fw is None!")
                fw.write(fragment)
                chrom_current = fragment.chrom
            if fw is not None:
                fw.close()
            f.close()

    @staticmethod
    def merge(args: argparse.Namespace) -> None:
        """Merge several FAM files into a single one.

        The last positional argument is the output file, the preceding
        ones are the inputs.
        """
        inputs = cast(list[str], args.inputs)
        if len(inputs) < 2:
            raise RuntimeError("At least one input file and one output file are required.")
        outfam = cast(str, args.outfam)

        # Reuse the header of the first input as the output template.
        fw = None
        for infam in inputs:
            f = FamFile(infam)
            if fw is None:
                fw = FamFile(outfam, "wb", f)
            for frag in f.fetch():
                fw.write(frag)
            f.close()
        if fw is not None:
            fw.close()
        

    # @staticmethod
    # def maskevent(args: argparse.Namespace) -> None:
    #     """Remove mismatch events overlapping the supplied SNP sites.

    #     The positions of the SNP file (BED/VCF) are used to blank the
    #     matching ``ME`` entries of each mate; the result is written to
    #     the ``CE`` tag.
    #     """
    #     with FamFile(args.in_fam) as f, FamFile(args.out_fam, "wb", f) as fw:
    #         if args.snps.endswith(".bed") or args.snps.endswith(".bed.gz"):
    #             snps = BedFile(args.snps)
    #         else:
    #             snps = VcfFile(args.snps)
    #         loader = ShiftLoader(snps)
    #         for fragment in f:
    #             positions = {snp.start for snp in loader.fetch(obj=fragment)}
    #             for align in (fragment.mate1, fragment.mate2):
    #                 segment = align.segment
    #                 if not segment.has_tag("ME"):
    #                     continue
    #                 events = []
    #                 for item in segment.get_tag("ME").split(";"):
    #                     if item == "":
    #                         continue
    #                     e = item.split(",")
    #                     if e[1] == "-" or e[2] == "-":
    #                         continue
    #                     if int(e[0]) in positions:
    #                         continue
    #                     events.append(item)
    #                 segment.set_tag("CE", ";".join(events))
    #             fw.write(fragment)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="A toolkit to process FAM file")
    subparsers = parser.add_subparsers(title="Available subcommands",
                                       dest="command")

    build_parser = subparsers.add_parser(
        "build", help="Build FAM file",
        description="Build a FAM file from a BAM file. Records must be pair-adjacent (e.g. query-name sorted).")
    _ = build_parser.add_argument("-i", "--input", help="Input BAM file")
    _ = build_parser.add_argument("-o", "--output", help="Output FAM file")

    split_parser = subparsers.add_parser(
        "split", help="Split FAM file by seqname",
        description="Split a FAM file into one FAM file per reference.")
    _ = split_parser.add_argument("-i", "--infam", help="Input FAM file")
    _ = split_parser.add_argument("-o", "--outdir", help="Output directory")
    _ = split_parser.add_argument(
        "-a", "--all", action="store_true",
        help="Output a file for every reference (even empty ones).")

    merge_parser = subparsers.add_parser(
        "merge", help="Merge splitted FAM files",
        description="Merge several FAM files into one.")
    _ = merge_parser.add_argument(
        "inputs", nargs="+",
        help="At least two input files")
    _ = merge_parser.add_argument("-o", "--output", help="Output FAM file")

    maskevent_parser = subparsers.add_parser(
        "maskevent", help="Mask events by SNP sites",
        description="Remove mismatch events overlapping the SNP sites (BED/VCF), writing the result to the CE tag.")
    _ = maskevent_parser.add_argument("in_fam", help="Input FAM file")
    _ = maskevent_parser.add_argument(
        "snps", help="SNP file (BED or VCF)")
    _ = maskevent_parser.add_argument("out_fam", help="Output FAM file")

    args = parser.parse_args()
    command = cast(str | None, args.command)
    if command is None:
        parser.print_help()
        return
    dispatch = {
        "build": FamTools.build,
        "split": FamTools.split,
        "merge": FamTools.merge,
        # "maskevent": FamTools.maskevent,
    }
    dispatch[command](args)


if __name__ == "__main__":
    main()
