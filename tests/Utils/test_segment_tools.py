#!/usr/scripts/env python

import importlib.util
import re
import tempfile
import unittest

import pysam

from pyBioInfo.Utils import SegmentTools


# Reverse flags used to build mate pairs and strand orientation.
FLAG_READ1 = 0x41   # 1 + 64: paired, read1
FLAG_READ2 = 0x81   # 1 + 128: paired, read2
FLAG_REVERSE = 0x10  # reverse strand


class TestSegmentTools(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.header = pysam.AlignmentHeader.from_dict({
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 2000}, {"SN": "chr2", "LN": 2000}],
        })

    def make_seg(self, name, cigar, ref_start=100, qseq=None, flag=0,
                 ref_id=0, mate_ref=None, mate_start=None):
        """Build a synthetic aligned segment attached to the shared header."""
        seg = pysam.AlignedSegment(self.header)
        seg.query_name = name
        seg.flag = flag
        seg.reference_id = ref_id
        seg.reference_start = ref_start
        seg.mapping_quality = 60
        seg.cigarstring = cigar
        qlen = sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
                   if op in "MIS=X")
        seg.query_sequence = "ACGT" * 50 if qseq is None else qseq
        seg.query_sequence = seg.query_sequence[:qlen].upper()
        seg.query_qualities = pysam.qualitystring_to_array("I" * qlen)
        if mate_ref is not None:
            seg.next_reference_id = mate_ref
            seg.next_reference_start = mate_start
        return seg

    # cmp / sort ------------------------------------------------------

    def test_cmp(self):
        """Ordering by chrom, then start, end and name."""
        a = self.make_seg("r2", "10M", 100)
        b = self.make_seg("r1", "10M", 50)
        c = self.make_seg("r3", "10M", 100)
        d = self.make_seg("r1", "10M", 100, ref_id=1)  # chr2
        self.assertEqual(SegmentTools.cmp(b, a), -1)
        self.assertEqual(SegmentTools.cmp(a, b), 1)
        self.assertEqual(SegmentTools.cmp(a, c), -1)  # same coords, name differs
        self.assertEqual(SegmentTools.cmp(a, a), 0)
        self.assertEqual(SegmentTools.cmp(d, a), 1)  # chr2 sorts after chr1

    def test_sort(self):
        """Coordinate sorting and name sorting."""
        a = self.make_seg("r2", "10M", 100)
        b = self.make_seg("r1", "10M", 50)
        c = self.make_seg("r3", "10M", 100)
        by_coord = SegmentTools.sort([a, b, c])
        self.assertEqual([s.query_name for s in by_coord], ["r1", "r2", "r3"])
        by_name = SegmentTools.sort([a, b, c], by_coordinates=False, by_name=True)
        self.assertEqual([s.query_name for s in by_name], ["r1", "r2", "r3"])
        self.assertRaises(RuntimeError, SegmentTools.sort, [a], by_coordinates=True, by_name=True)
        self.assertRaises(RuntimeError, SegmentTools.sort, [a], by_coordinates=False, by_name=False)

    # blocks / lengths ------------------------------------------------

    def test_get_blocks(self):
        """Block extraction handles N skips and D operations."""
        s = self.make_seg("r1", "5M3N10M", 100)
        self.assertEqual(SegmentTools.get_blocks(s), [[100, 105], [108, 118]])
        s2 = self.make_seg("r2", "5M2D10M", 100)
        self.assertEqual(SegmentTools.get_blocks(s2), [[100, 117]])
        self.assertEqual(SegmentTools.get_blocks(s2, False), [[100, 105], [107, 117]])
        # a deletion at the very start is filled and merged with the M block.
        s3 = self.make_seg("r3", "2D5M", 100)
        self.assertEqual(SegmentTools.get_blocks(s3), [[100, 107]])
        self.assertEqual(SegmentTools.get_blocks(s3, False), [[102, 107]])

    def test_get_mapped_length(self):
        """Mapped length counts M and (optionally) D operations."""
        s = self.make_seg("r1", "5M2D10M", 100)
        self.assertEqual(SegmentTools.get_mapped_length(s), 17)
        self.assertEqual(SegmentTools.get_mapped_length(s, False), 15)

    def test_get_aligned_sequence(self):
        """Deletions become '-'/-1, insertions and clips are dropped."""
        s = self.make_seg("r1", "3S4M2D3M1I2M", 100)
        seq, quals = SegmentTools.get_aligned_sequence(s)
        self.assertEqual(len(seq), len(quals))
        self.assertEqual(seq, "TACG--TACTA")  # 4M, 2D, 3M, (1I dropped), 2M
        self.assertEqual(quals[4], -1)
        self.assertEqual(quals[5], -1)
        # Insertion bases do not consume the aligned reference length.
        for q in quals:
            if q != -1:
                self.assertEqual(q, 40)  # “I” quality

    # parsed CIGAR ----------------------------------------------------

    def test_parse_cigar_from_segment(self):
        """Parse a segment CIGAR via parse_cigar(segment=...)."""
        s = self.make_seg("r1", "5M2D10M", 100)
        parsed = SegmentTools.parse_cigar(s)
        self.assertEqual(parsed, [
            ("M", (0, 5), (0, 5), (100, 105)),
            ("D", (5, 7), (5, 5), (105, 107)),
            ("M", (7, 17), (5, 15), (107, 117)),
        ])

    @unittest.skipUnless(importlib.util.find_spec("cigar"),
                         "optional 'cigar' library is not installed")
    def test_parse_cigar_from_cigarstring(self):
        """Parse a CIGAR string via parse_cigar(start=..., cigarstring=...)."""
        parsed = SegmentTools.parse_cigar(start=100, cigarstring="5M2D10M")
        self.assertEqual(parsed, [
            ("M", (0, 5), (0, 5), (100, 105)),
            ("D", (5, 7), (5, 5), (105, 107)),
            ("M", (7, 17), (5, 15), (107, 117)),
        ])
        # Hard clip consumes the read indices but not the reference.
        parsed2 = SegmentTools.parse_cigar(start=100, cigarstring="3H5M2I4M")
        self.assertEqual(parsed2[0], ("H", (0, 0), (0, 3), (100, 100)))
        self.assertEqual(parsed2[1], ("M", (0, 5), (3, 8), (100, 105)))

    @unittest.skipUnless(importlib.util.find_spec("cigar"),
                         "optional 'cigar' library is not installed")
    def test_parse_cigar_paths_consistent(self):
        """The segment and string entry points must agree (mapped indices
        accumulate correctly across D/N operations)."""
        for cigar in ("5M2D10M", "3S4M5N2M", "3H5M2I4M"):
            s = self.make_seg("r1", cigar, 100)
            from_seg = SegmentTools.parse_cigar(s)
            from_str = SegmentTools.parse_cigar(start=100, cigarstring=cigar)
            self.assertEqual(from_seg, from_str)

    def test_parse_cigar_dispatch(self):
        """parse_cigar validates its input arguments."""
        s = self.make_seg("r1", "5M", 100)
        self.assertRaises(RuntimeError, SegmentTools.parse_cigar)
        self.assertRaises(RuntimeError, SegmentTools.parse_cigar,
                          segment=s, cigarstring="5M")
        # A cigar string without a start position is invalid.
        self.assertRaises(RuntimeError, SegmentTools.parse_cigar,
                          cigarstring="5M")

    # MD tag ----------------------------------------------------------

    def test_parse_md_tag(self):
        """Parse MD tags into match/mismatch/deletion entries."""
        s = self.make_seg("r1", "10M", 100)
        s.set_tag("MD", "3A6")
        self.assertEqual(SegmentTools.parse_md_tag(s), [
            ("=", (0, 3), None),
            ("X", (3, 4), "A"),
            ("=", (4, 10), None),
        ])
        s.set_tag("MD", "5^ACG4")
        self.assertEqual(SegmentTools.parse_md_tag(s), [
            ("=", (0, 5), None),
            ("D", (5, 8), "ACG"),
            ("=", (8, 12), None),
        ])

    # clipped bases ---------------------------------------------------

    def test_get_clipped(self):
        """Clipped bases at both ends, honouring strand orientation."""
        s1 = self.make_seg("r1", "3S4M", 100)
        self.assertEqual(SegmentTools.get_clipped(s1), (3, 0))
        self.assertEqual(SegmentTools.get_clipped(s1, soft_clip=False), (0, 0))
        s2 = self.make_seg("r2", "2H3S4M1H", 100)
        self.assertEqual(SegmentTools.get_clipped(s2), (2, 1))
        self.assertEqual(SegmentTools.get_clipped(s2, hard_clip=False), (0, 0))
        s3 = self.make_seg("r3", "3S4M", 100, flag=FLAG_REVERSE)
        self.assertEqual(SegmentTools.get_clipped(s3), (3, 0))
        self.assertEqual(SegmentTools.get_clipped(s3, stranded=True), (0, 3))

    # query base ------------------------------------------------------

    def test_get_query_base(self):
        """Read base aligned at a reference position."""
        s = self.make_seg("r1", "5M2D10M", 100, qseq="AAAAACCCCCGGGGG")
        self.assertEqual(SegmentTools.get_query_base(s, 100), "A")
        self.assertEqual(SegmentTools.get_query_base(s, 104), "A")
        self.assertEqual(SegmentTools.get_query_base(s, 105), "-")  # deletion
        self.assertEqual(SegmentTools.get_query_base(s, 106), "-")  # deletion
        self.assertEqual(SegmentTools.get_query_base(s, 107), "C")
        self.assertEqual(SegmentTools.get_query_base(s, 108), "C")
        # reference skip (N) yields None; outside range yields None.
        s2 = self.make_seg("r2", "5M3N5M", 100, qseq="AAAACCCCCCC")
        self.assertIsNone(SegmentTools.get_query_base(s2, 105))
        self.assertIsNone(SegmentTools.get_query_base(s2, 90))

    # reference sequence ----------------------------------------------

    def test_get_reference_sequence(self):
        """Reference sequence rebuilt from query + MD tag."""
        s = self.make_seg("r1", "5M", 100, qseq="ACGTA")
        s.set_tag("MD", "5")
        self.assertEqual(SegmentTools.get_reference_sequence(s), "ACGTA")

    def test_get_reference_sequence_base(self):
        """Reference base occurrence counting (upper case)."""
        s = self.make_seg("r1", "5M", 100, qseq="ACGTA")
        s.set_tag("MD", "5")
        counts = SegmentTools.get_reference_sequence_base(s)
        self.assertEqual(dict(counts), {"A": 2, "C": 1, "G": 1, "T": 1})

    # grouping --------------------------------------------------------

    def test_load_same_name_segments_from_iterable(self):
        """Group a sequential stream by query name."""
        pool = [
            self.make_seg("g1", "10M", 100),
            self.make_seg("g1", "10M", 200),
            self.make_seg("g2", "10M", 300),
        ]
        groups = list(SegmentTools.load_same_name_segments(pool))
        self.assertEqual(len(groups), 2)
        self.assertEqual([g.query_name for g in groups[0]], ["g1", "g1"])
        self.assertEqual([g.query_name for g in groups[1]], ["g2"])

    def test_load_same_name_segments_from_file(self):
        """Group segments read from a SAM file."""
        sam = "\n".join([
            "@HD\tVN:1.6\tSO:queryname",
            "@SQ\tSN:chr1\tLN:2000",
            "g1\t0\tchr1\t101\t60\t5M\t*\t0\t0\tAAAAA\tIIIII",
            "g1\t0\tchr1\t201\t60\t5M\t*\t0\t0\tAAAAA\tIIIII",
            "g2\t0\tchr1\t301\t60\t5M\t*\t0\t0\tAAAAA\tIIIII",
        ])
        with tempfile.NamedTemporaryFile(mode="w", suffix=".sam",
                                         delete=False) as fh:
            fh.write(sam + "\n")
            path = fh.name
        try:
            groups = list(SegmentTools.load_same_name_segments(path))
            self.assertEqual(len(groups), 2)
            self.assertEqual([g.query_name for g in groups[0]], ["g1", "g1"])
            self.assertEqual([g.query_name for g in groups[1]], ["g2"])
        finally:
            import os
            os.remove(path)

    # mates -----------------------------------------------------------

    def test_is_mate(self):
        """Mate detection across name, positions, flags and MC tags."""
        s1 = self.make_seg("pair", "10M", 100, flag=FLAG_READ1,
                           mate_ref=0, mate_start=200)
        s2 = self.make_seg("pair", "5M", 200, flag=FLAG_READ2,
                           mate_ref=0, mate_start=100)
        s1.cigarstring = "10M"
        s2.cigarstring = "5M"
        s1.set_tag("MC", "5M")
        s2.set_tag("MC", "10M")
        self.assertTrue(SegmentTools.is_mate(s1, s2))
        # Different query names are never mates.
        s3 = self.make_seg("other", "5M", 200, flag=FLAG_READ2,
                           mate_ref=0, mate_start=100)
        self.assertFalse(SegmentTools.is_mate(s1, s3))
        # Both ends keep the same read flag (both read1) is rejected.
        s4 = self.make_seg("pair", "5M", 200, flag=FLAG_READ1,
                           mate_ref=0, mate_start=100)
        self.assertFalse(SegmentTools.is_mate(s1, s4))
        # One end secondary (flag 0x100) is rejected.
        s5 = self.make_seg("pair", "5M", 200, flag=FLAG_READ2 | 0x100,
                           mate_ref=0, mate_start=100)
        self.assertFalse(SegmentTools.is_mate(s1, s5))


if __name__ == "__main__":
    _ = unittest.main()
