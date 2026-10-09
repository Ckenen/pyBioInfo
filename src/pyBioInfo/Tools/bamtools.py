import argparse
import os
import shutil
import subprocess
from typing import cast
from pyBioInfo.Range import GRange
from pyBioInfo.Tools.bedtools import BedTools
from pyBioInfo.IO.File import BamFile, FamFile
import pysam
from collections import defaultdict
import pandas as pd
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class Add_NH_HI_Tag_Worker(object):
    @classmethod  
    def run_cmd(cls, args: list[str]) -> None:
        cmd = " ".join(args)
        logging.info("RUN CMD: %s" % cmd)
        _ = subprocess.run(args, check=True)

    @classmethod  
    def load_same_name_segments(cls, f: pysam.AlignmentFile):
        segments = None
        for segment in f.fetch(until_eof=True):
            if segments is None:
                segments = [segment]
            elif segment.query_name == segments[0].query_name:
                segments.append(segment)
            else:
                yield segments
                segments = [segment]
        if segments is not None:
            yield segments

    @classmethod  
    def detect_layout(cls, bam_file: str) -> str:
        with pysam.AlignmentFile(bam_file) as f:
            for s in f:
                if s.is_paired:
                    return "PE"
                else:
                    return "SE"
        return "SE"

    @classmethod  
    def output_summary(cls, counter: defaultdict[int, int], outfile: str) -> None:
        total = sum(counter.values())
        with open(outfile, "w+")  as fw:
            _ = fw.write("NH\tCount\tRatio\n")
            for k, v in sorted(counter.items()):
                _ = fw.write("%d\t%d\t%.6f\n" % (k, v, v / total))
            
    @classmethod  
    def add_nh_hi_tag_for_se(cls, input_bam: str, output_bam: str, summary_tsv: str) -> None:
        num_unmap = 0
        num_mapped = 0
        nh_counter: defaultdict[int, int] = defaultdict(int)
        with pysam.AlignmentFile(input_bam) as f, pysam.AlignmentFile(output_bam, "wb", template=f) as fw:
            for segments in cls.load_same_name_segments(f):
                segments_unmap: list[pysam.AlignedSegment] = []
                segments_mapped: list[pysam.AlignedSegment] = []
                for s in segments:
                    if s.is_unmapped:
                        segments_unmap.append(s)
                    else:
                        segments_mapped.append(s)
                num_unmap += len(segments_unmap)
                num_mapped += len(segments_mapped)
                            
                primary = None
                secondaries: list[pysam.AlignedSegment] = []
                for s in segments_mapped:
                    if s.is_secondary:
                        secondaries.append(s)
                    elif s.is_supplementary:
                        continue
                        # raise NotImplementedError()
                    else:
                        primary = s

                nh = len(segments)
                nh_counter[nh] += 1
                hi = 1
                if primary is not None:
                    primary.set_tag("NH", nh)
                    primary.set_tag("HI", hi)
                    hi += 1
                    _ = fw.write(primary)
                for s in sorted(secondaries, key=lambda x: [-x.mapping_quality, x.reference_name, x.reference_start, x.reference_end]):
                    s.set_tag("NH", nh)
                    s.set_tag("HI", hi)
                    hi += 1
                    _ = fw.write(s)
        
        cls.output_summary(nh_counter, summary_tsv)
        
        num_total = num_unmap + num_mapped
        logging.info("Number of segments: %d" % num_total)
        logging.info("Number of unmapped segments: %d (%.6f)" % (num_unmap, num_unmap / num_total))
        logging.info("Number of mapped segments: %d (%.6f)" % (num_mapped, num_mapped / num_total))

    @classmethod  
    def add_nh_hi_tag_for_pe(cls, input_bam: str, output_bam: str, proper_pair: bool, summary_tsv: str) -> None:
        num_unmap = 0
        num_mate_unmap = 0
        num_both_mapped = 0
        num_proper_mapped = 0
        num_wrote_pairs = 0
        nh_counter: defaultdict[int, int] = defaultdict(int)
        with pysam.AlignmentFile(input_bam) as f, pysam.AlignmentFile(output_bam, "wb", template=f) as fw:
            for segments in cls.load_same_name_segments(f):
                # if i >= 1000:
                #     break
                read_name = segments[0].query_name
                segments_unmap: list[pysam.AlignedSegment] = []
                segments_mate_unmap: list[pysam.AlignedSegment] = []
                segments_both_mapped: list[pysam.AlignedSegment] = []
                segments_proper_mapped: list[pysam.AlignedSegment] = []
                for s in segments:
                    if s.is_unmapped:
                        segments_unmap.append(s)
                    elif s.mate_is_unmapped:
                        segments_mate_unmap.append(s)
                    else:
                        segments_both_mapped.append(s)
                        if s.is_proper_pair:
                            segments_proper_mapped.append(s)
                            
                num_unmap += len(segments_unmap)
                num_mate_unmap += len(segments_mate_unmap)
                num_both_mapped += len(segments_both_mapped)
                num_proper_mapped += len(segments_proper_mapped)
                
                if proper_pair:
                    segments_both_mapped = segments_proper_mapped
                    
                data1 = defaultdict(list)
                data2 = defaultdict(list)
                for s in segments_both_mapped:
                    chrom1 = s.reference_name
                    start1 = s.reference_start
                    chrom2 = s.next_reference_name
                    start2 = s.next_reference_start
                    if chrom1 is None or chrom2 is None:
                        raise ValueError(f"chrom1 is None or chrom2 is None: {chrom1}, {chrom2}")
                        
                    is_primary = not (s.is_secondary or s.is_supplementary) # is primary
                    if s.has_tag("MC"):
                        cigar1 = s.cigarstring
                        if not isinstance(cigar1, str):
                            raise TypeError(f"cigar1 is not str: {cigar1}")
                        cigar2 = cast(str, s.get_tag("MC"))

                        if s.is_read1:
                            data1[(chrom1, start1, chrom2, start2, is_primary, cigar1, cigar2)].append(s)
                        elif s.is_read2:
                            data2[(chrom2, start2, chrom1, start1, is_primary, cigar2, cigar1)].append(s)
                        else:
                            assert False
                    else:
                        if s.is_read1:
                            data1[(chrom1, start1, chrom2, start2, is_primary)].append(s)
                        elif s.is_read2:
                            data2[(chrom2, start2, chrom1, start1, is_primary)].append(s)
                        else:
                            assert False
                        
                for k, v in data1.items():
                    if len(v) > 1:
                        logging.warning("Warning: one key map to multiple segments (Mate: 1, Name: %s, Size: %d, key: %s)" % (read_name, len(v), " ".join(map(str, k))))
                for k, v in data2.items():
                    if len(v) > 1:
                        logging.warning("Warning: one key map to multiple segments (Mate: 2, Name: %s, Size: %d, key: %s)" % (read_name, len(v), " ".join(map(str, k))))
                        
                ks1, ks2 = set(data1.keys()), set(data2.keys())
                ks_common = ks1 & ks2
                
                if len(ks_common) == 0:
                    continue
                
                # TODO
                # ks1_only = ks1 - ks2
                # ks2_only = ks2 - ks1
                # print(len(ks_common), len(ks1_only), len(ks2_only))
                # assert len(ks1_only) == 0
                # assert len(ks2_only) == 0
                    
                array = []
                for k in ks_common:
                    # assert len(data1[k]) == 1 and len(data2[k]) == 1
                    s1, s2 = data1[k][0], data2[k][0]
                    array.append([[-s1.mapping_quality, -s2.mapping_quality, k], [s1, s2]])
                nh = len(array)
                nh_counter[nh] += 1
                primary, secondaries = None, []
                for item in array:
                    # k = item[0][2]
                    # is_primary = k[4]
                    s1, s2 = item[1]
                    if s1.is_secondary and s2.is_secondary:
                        secondaries.append(item)
                    # elif not s1.is_secondary and not s2.is_secondary and not s1.is_supplementary and not s2.is_supplementary:
                    elif item[0][2][4]: # is_primary
                        primary = item
                    else:
                        assert False
                hi = 1
                if primary is not None:
                    s1, s2 = primary[1]
                    s1.set_tag("NH", nh)
                    s1.set_tag("HI", hi)
                    s2.set_tag("NH", nh)
                    s2.set_tag("HI", hi)
                    hi += 1
                    num_wrote_pairs += 1
                    _ = fw.write(s1)
                    _ = fw.write(s2)
                for item in sorted(secondaries, key=lambda x: x[0]):
                    s1, s2 = item[1]
                    s1.set_tag("NH", nh)
                    s1.set_tag("HI", hi)
                    s2.set_tag("NH", nh)
                    s2.set_tag("HI", hi)
                    hi += 1
                    num_wrote_pairs += 1
                    _ = fw.write(s1)
                    _ = fw.write(s2)
        
        cls.output_summary(nh_counter, summary_tsv)
        
        num_total = num_unmap + num_mate_unmap + num_both_mapped
        logging.info("Number of segments: %d" % num_total)
        logging.info("Number of unmapped segments: %d (%.6f)" % (num_unmap, num_unmap / num_total))
        logging.info("Number of mate unmapped segments: %d (%.6f)" % (num_mate_unmap, num_mate_unmap / num_total))
        logging.info("Number of both mapped segments: %d (%.6f)" % (num_both_mapped, num_both_mapped / num_total))
        logging.info("Number of proper-paired segments: %d (%.6f)" % (num_proper_mapped, num_proper_mapped / num_total))
        logging.info("Number of output paired segments: %d (%.6f)" % (num_wrote_pairs * 2, num_wrote_pairs * 2 / num_total))
        
    @classmethod            
    def add_nh_hi_tag(cls, input_bam: str, output_bam: str, proper_pair: bool, summary_tsv: str) -> None:
        layout = cls.detect_layout(input_bam)
        if layout == "SE":
            cls.add_nh_hi_tag_for_se(input_bam, output_bam, summary_tsv)
        elif layout == "PE":
            cls.add_nh_hi_tag_for_pe(input_bam, output_bam, proper_pair, summary_tsv)
        else:
            logging.error("Unknown layout: %s" % layout)
            exit(1)
    
    
class BamTools(object):
    @classmethod
    def bam_to_bed(cls, args: argparse.Namespace): 
        bamfile = cast(str, args.input)
        bedfile = cast(str, args.output)
        bedfile2 = bedfile + ".gz"
        merge_frag = cast(bool, args.merge_fragment)
        fill_frag_gap = cast(bool, args.fill_fragment_gap)
        build_index = cast(bool, args.build_index)
        if fill_frag_gap:
            merge_frag = True
            # assert merge_frag, "--fill-fragment-gap requires --merge-fragment"
        
        logging.info(f"Input BAM: {bamfile}")
        logging.info(f"Output BED: {bedfile}")
        logging.info(f"Output BED.gz: {bedfile2}")
        logging.info(f"Merge fragment: {merge_frag}")
        logging.info(f"Fill fragment gap: {fill_frag_gap}")
        logging.info(f"Build index: {build_index}")
        
        layout = BamFile.detect_layout(bamfile)

        with open(bedfile, "w+") as fw:
            if layout =="SE":
                assert bamfile.endswith(".bam")
                with BamFile(bamfile, random=True) as f:
                    for align in f.fetch():
                        segment = align.segment
                        if segment.has_tag("NH") and segment.has_tag("HI"):
                            align.name = "%s[NH=%d][HI=%d]" % (align.name, segment.get_tag("NH"), segment.get_tag("HI"))
                        _ = fw.write(align.format("BED") + "\n")
            else:
                if merge_frag:
                    assert bamfile.endswith(".bam") or bamfile.endswith(".fam")
                    rd = bamfile.endswith(".bam")
                    if fill_frag_gap:
                        with FamFile(bamfile, random=rd) as f:
                            for fragment in f.fetch():
                                mate1 = fragment.mate1
                                segment1 = mate1.segment
                                mate2 = fragment.mate2

                                if segment1.has_tag("NH") and segment1.has_tag("HI"):
                                    fragment.name = "%s[NH=%d][HI=%d]" % (fragment.name, segment1.get_tag("NH"), segment1.get_tag("HI"))

                                if mate1.is_overlap(mate2):
                                    _ = fw.write(fragment.format("BED") + "\n")
                                else:
                                    blocks1 = [[x, y] for x, y in mate1.blocks]
                                    blocks2 = [[x, y] for x, y in mate2.blocks]
                                    if blocks1[-1][1] <= blocks2[0][0]:
                                        pass
                                    elif blocks2[-1][1] <= blocks1[0][0]:
                                        blocks1, blocks2 = blocks2, blocks1
                                    else:
                                        raise RuntimeError("Cannot fill gap between paired alignments")
                                    blocks1[-1][1] = blocks2[0][1]
                                    blocks = blocks1 + blocks2[1:]
                                    obj = GRange(chrom=fragment.chrom, name=fragment.name, strand=fragment.strand, blocks=blocks)
                                    _ = fw.write(obj.format("BED") + "\n")
                    else:
                        with FamFile(bamfile, random=rd) as f:
                            for fragment in f.fetch():
                                mate1 = fragment.mate1
                                segment1 = mate1.segment
                                if segment1.has_tag("NH") and segment1.has_tag("HI"):
                                    fragment.name = "%s[NH=%d][HI=%d]" % (fragment.name, segment1.get_tag("NH"), segment1.get_tag("HI"))
                                _ = fw.write(fragment.format("BED") + "\n")
                else:
                    assert bamfile.endswith(".bam")
                    with BamFile(bamfile, random=True) as f:
                        for align in f.fetch():
                            segment = align.segment
                            if segment.has_tag("NH") and segment.has_tag("HI"):
                                align.name = "%s[NH=%d][HI=%d]" % (align.name, segment.get_tag("NH"), segment.get_tag("HI"))
                            align.name = "%s[RI=%s]" % (align.name, "1" if align.segment.is_read1 else "2")
                            _ = fw.write(align.format("BED") + "\n")

        if build_index:
            _ = subprocess.run(["bgzip", "-f", bedfile, "-o", bedfile2], check=True)
            _ = subprocess.run(["tabix", "-p", "bed", bedfile2], check=True)

        logging.info("Done")

    @classmethod
    def bam_to_bigwig(cls, args: argparse.Namespace) -> None:
        bamfile = cast(str, args.bamfile)
        outfile1 = cast(str, args.outfile1)
        outfile2 = cast(str, args.outfile2)
        scale = cast(float, args.scale)
        stranded = cast(bool, args.stranded)
        merge_fragment = cast(bool, args.merge_frag)
        fill_fragment_gap: float = cast(float, args.fill_frag_gap)

        logging.info(f"Input BAM: {bamfile}")
        logging.info(f"Output BigWig: {outfile1}")
        logging.info(f"Output BigWig (reverse strand): {outfile2}")
        logging.info(f"Scale: {scale}")
        logging.info(f"Stranded: {stranded}")
        logging.info(f"Merge fragment: {merge_fragment}")
        logging.info(f"Fill fragment gap: {fill_fragment_gap}")

        tmpdir = outfile1 + ".bamtools_bam2bw_tmp"
        os.makedirs(tmpdir)
        bedfile = os.path.join(tmpdir, "converted.bed")
        _ = cls.bam_to_bed(args=argparse.Namespace(input=bedfile, output=bedfile, merge_frament=merge_fragment, fill_fragment_gap=fill_fragment_gap))
        genomesize = os.path.join(tmpdir, "genome_sizes.txt")
        _ = cls.extract_genome_sizes(args=argparse.Namespace(input=bamfile, output=genomesize))
        _ = BedTools.bed_to_bigwig(args=argparse.Namespace(input=bedfile, output1=outfile1, output2=outfile2, genome_sizes=genomesize, scale=scale, stranded=stranded))
        shutil.rmtree(tmpdir)

        logging.info("Done")
        
    @classmethod
    def extract_genome_sizes(cls, args: argparse.Namespace) -> None:
        bamfile = cast(str, args.input)
        outfile = cast(str, args.output)
        with BamFile(bamfile, random=False) as f, open(outfile, "w+") as fw:
            references = f.references
            if references is None:
                raise RuntimeError("Failed to get reference information from BAM file")
            for chrom, length in references.items():
                _ = fw.write("%s\t%d\n" % (chrom, length))

    @classmethod
    def add_nh_hi_tag(cls, args: argparse.Namespace):
        bamfile = cast(str, args.input)
        outfile = cast(str, args.output)
        threads = cast(int, args.threads)
        input_is_sorted = cast(bool, args.input_is_sorted)
        sort_output = cast(bool, args.sort_output)
        proper_pair = cast(bool, args.proper_pair)

        logging.info(f"Input BAM: {bamfile}")
        logging.info(f"Output BAM: {outfile}")
        logging.info(f"Threads: {threads}")
        logging.info(f"Input is sorted: {input_is_sorted}")
        logging.info(f"Sort output: {sort_output}")
        logging.info(f"Proper pair: {proper_pair}")

        tmpdir = outfile + ".bamtools_add_nh_hi_tag_tmp"

        Add_NH_HI_Tag_Worker.run_cmd(["mkdir", "-p", tmpdir])
        # os.makedirs(tmpdir, exist_ok=True)

        bamfile1 = bamfile
        bamfile2 = bamfile
        bamfile3 = outfile
        bamfile4 = outfile
        summary_tsv = outfile + ".tsv"

        if not input_is_sorted:
            bamfile2 = os.path.join(tmpdir, "input.sorted_by_name.bam")

        if sort_output:
            bamfile3 = os.path.join(tmpdir, "output.sorted_by_name.bam")

        if bamfile1 != bamfile2:
            Add_NH_HI_Tag_Worker.run_cmd(["samtools", "sort", "-@", str(threads), "-N", "-o", bamfile2, bamfile1])

        Add_NH_HI_Tag_Worker.add_nh_hi_tag(bamfile2, bamfile3, proper_pair, summary_tsv)

        if bamfile3 != bamfile4:
            Add_NH_HI_Tag_Worker.run_cmd(["samtools", "sort", "-@", str(threads), "-o", bamfile4, bamfile3])
            Add_NH_HI_Tag_Worker.run_cmd(["samtools", "index", "-@", str(threads), bamfile4])

        Add_NH_HI_Tag_Worker.run_cmd(["rm", "-rf", tmpdir])

        logging.info("Done")

    @classmethod
    def stat_chrom_reads(cls, args: argparse.Namespace) -> None:
        bamfile = cast(str, args.input)
        outfile = cast(str, args.output)
        primary_only = cast(bool, args.primary_only)
        proper_pair_only = cast(bool, args.proper_pair_only)

        logging.info(f"Input BAM: {bamfile}")
        logging.info(f"Output TSV: {outfile}")
        logging.info(f"Require primary only: {primary_only}")
        logging.info(f"Require proper pair: {proper_pair_only}")

        rows = []

        with pysam.AlignmentFile(bamfile) as f:

            for chrom in f.references:

                n1, n2 = 0, 0

                for s in f.fetch(chrom):
                    if primary_only:
                        if s.is_secondary:
                            continue
                        if s.is_supplementary:
                            continue

                    if s.is_paired:
                        if proper_pair_only and not s.is_proper_pair:
                            continue
                        if s.is_read1:
                            if s.is_reverse:
                                n2 += 1
                            else:
                                n1 += 1
                        elif s.is_read2:
                            if s.is_reverse:
                                n1 += 1
                            else:
                                n2 += 1
                        else:
                            assert False
                    else:
                        if s.is_reverse:
                            n2 += 1
                        else:
                            n1 += 1

                rows.append([chrom, n1, n2, n1 + n2])

        df = pd.DataFrame(rows)
        df.columns = ["Chrom", "Forward", "Reverse", "Both"]
        df["Total"] = df["Both"].sum()
        df["Proportion"] = df["Both"] / df["Total"]
        df["File"] = bamfile
        df.to_csv(outfile, sep="\t" if outfile.endswith(".tsv") else ",", index=False)

        logging.info("Done")

    # @classmethod
    # def extract(cls, args):
    #     usage = "bamtools.py extract [options] input.bam output.bam"
    #     parser = OptionParser(usage=usage)
    #     parser.add_option("-c", "--chrom", dest="chrom", type="string")
    #     parser.add_option("-s", "--start", dest="start",
    #                       default=None, type="int")
    #     parser.add_option("-e", "--end", dest="end", default=None, type="int")
    #     options, args = parser.parse_args(args)

    #     infile, outfile = args
    #     chrom = options.chrom
    #     start = options.start
    #     end = options.end
    #     assert chrom is not None

    #     with pysam.AlignmentFile(infile) as f, pysam.AlignmentFile(outfile, "wb", template=f) as fw:
    #         for segment in f.fetch(contig=chrom, start=start, end=end):
    #             fw.write(segment)

    # @classmethod
    # def getevent(cls, args):
    #     usage = "bamtools.py getevent input.bam output.bam"
    #     parser = OptionParser(usage=usage)
    #     options, args = parser.parse_args(args)
    #     infile, outfile = args
    #     with pysam.AlignmentFile(infile) as f, pysam.AlignmentFile(outfile, "wb", f) as fw:
    #         for segment in f:
    #             if segment.is_unmapped or segment.is_secondary or segment.is_supplementary:
    #                 continue
    #             bases = SegmentTools.get_reference_sequence_base(segment)
    #             rc = ";".join(["%s,%d" % (k, bases[k]) for k in sorted(bases.keys())])
    #             items = []
    #             for e in SegmentTools.get_events(segment):
    #                 if isinstance(e[3], list):
    #                     e[3] = "/".join(map(str, e[3]))
    #                 items.append(",".join(map(str, e)))
    #             me = ";".join(items)
    #             segment.set_tag("RC", rc)
    #             segment.set_tag("ME", me)
    #             fw.write(segment)
                
    # @classmethod
    # def maskevent(cls, args):
    #     usage = "bamtools.py maskevent input.bam snps.bed output.bam"
    #     parser = OptionParser(usage=usage)
    #     options, args = parser.parse_args(args)
    #     infile1, infile2, outfile = args
        
    #     with BamFile(infile1) as f, BamFile(outfile, "wb", f) as fw:
    #         if infile2.endswith(".bed") or infile2.endswith(".bed.gz"):
    #             snps = BedFile(infile2)
    #         elif infile2.endswith(".vcf") or infile2.endswith(".vcf.gz"):
    #             snps = VcfFile(infile2)
    #         loader = ShiftLoader(snps)
    #         for align in f:
    #             positions = set([snp.start for snp in loader.fetch(obj=align)])
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
    #             fw.write(align)

