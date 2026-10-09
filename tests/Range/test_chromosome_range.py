# Testing module of pyBioInfo.Range.chromosome_range

import unittest

from pyBioInfo.Range import CRange, IRange


class TestCRange(unittest.TestCase):
    """Testing class of pyBioInfo.Range.chromosome_range.CRange."""

    def test_constructor(self):
        """Construction with the chrom, start and end attributes."""
        obj = CRange("chr1", 5, 15)
        self.assertEqual("chr1", obj.chrom)
        self.assertEqual(5, obj.start)
        self.assertEqual(15, obj.end)
        # Invalid coordinates are rejected (inherited from IRange).
        self.assertRaises(ValueError, CRange, "chr1", -1, 10)
        self.assertRaises(ValueError, CRange, "chr1", 20, 10)

    def test_chrom_is_plain_attribute(self):
        """The chrom attribute is writable, unlike start and end."""
        obj = CRange("chr1", 5, 15)
        obj.chrom = "chr2"
        self.assertEqual("chr2", obj.chrom)
        self.assertRaises(AttributeError, setattr, obj, "start", 8)
        self.assertRaises(AttributeError, setattr, obj, "end", 18)

    def test_inherited_protocol(self):
        """len(), str() and the sequence protocol are inherited from IRange."""
        obj = CRange("chr1", 10, 20)
        self.assertEqual(10, len(obj))
        self.assertEqual("CRange: 10-20", str(obj))
        self.assertEqual([10, 20], list(iter(obj)))
        self.assertEqual(10, obj[0])
        self.assertEqual(20, obj[1])

    def test_strand(self):
        """Strand accessor, validation and direction conveniences."""
        # Valid strands.
        for strand in (None, ".", "+", "-"):
            self.assertIsInstance(CRange("chr1", 10, 20, strand=strand), CRange)
        # Invalid strand (validated by the setter).
        self.assertRaises(ValueError, CRange, "chr1", 10, 20, strand="S")
        obj = CRange("chr1", 10, 20)
        self.assertRaises(ValueError, setattr, obj, "strand", "S")
        # Direction conveniences.
        obj = CRange("chr1", 10, 20, strand="+")
        self.assertEqual("+", obj.strand)
        self.assertTrue(obj.is_forward)
        self.assertFalse(obj.is_reverse)
        obj.strand = "-"
        self.assertEqual("-", obj.strand)
        self.assertFalse(obj.is_forward)
        self.assertTrue(obj.is_reverse)
        obj.strand = None
        self.assertIsNone(obj.strand)
        self.assertFalse(obj.is_forward)
        self.assertFalse(obj.is_reverse)

    def test_comparison(self):
        """Ordering groups by chromosome first, then by coordinates."""
        obj1 = CRange("chr1", 10, 20)
        obj2 = CRange("chr1", 15, 25)
        obj3 = CRange("chr2", 5, 15)
        self.assertLess(obj1, obj2)
        self.assertLess(obj2, obj3)
        self.assertGreater(obj3, obj1)
        # Same chrom and start: ordered by the end coordinate.
        self.assertLess(CRange("chr1", 10, 20), CRange("chr1", 10, 30))
        # Equality of the sort key.
        self.assertEqual(obj1, CRange("chr1", 10, 20))
        self.assertNotEqual(obj1, obj2)
        self.assertNotEqual(obj1, 123)
        # Plain IRange operands fall back to the coordinate-only comparison.
        self.assertGreater(obj1, IRange(5, 15))
        self.assertLess(obj1, IRange(15, 25))
        self.assertEqual(obj1, IRange(10, 20))

    def test_index_and_position_strandness(self):
        """Strand-aware index/position conversions."""
        # Forward strand: plain coordinate semantics.
        obj = CRange("chr1", 10, 20, strand="+")
        self.assertEqual(2, obj.get_index(12))
        self.assertEqual(12, obj.get_position(2))
        # Reverse strand: index runs from the end towards the start.
        obj = CRange("chr1", 10, 20, strand="-")
        self.assertEqual(0, obj.get_index(19))   # end - 1 is the first base
        self.assertEqual(9, obj.get_index(10))
        self.assertEqual(19, obj.get_position(0))
        self.assertEqual(10, obj.get_position(9))
        self.assertEqual(10, obj.get_position(-1))  # negative counts from the strand-aware end
        self.assertRaises(ValueError, obj.get_position, 10)
        self.assertRaises(ValueError, obj.get_position, -11)
        # strandness=False disables the reversal.
        self.assertEqual(2, obj.get_index(12, strandness=False))
        self.assertEqual(12, obj.get_position(2, strandness=False))
        # No strand / unknown strand: plain semantics.
        for strand in (None, "."):
            obj = CRange("chr1", 10, 20, strand=strand)
            self.assertEqual(2, obj.get_index(12))
            self.assertEqual(12, obj.get_position(2))

    def test_distance_relationships(self):
        """Cross-chromosome ranges are always away; same-chrom uses coordinates."""
        obj = CRange("chr1", 10, 20)
        # Same chrom: normal coordinate-based relations.
        self.assertTrue(obj.is_overlap(CRange("chr1", 15, 25)))
        self.assertTrue(obj.is_adjacent(CRange("chr1", 20, 25)))
        self.assertTrue(obj.is_away(CRange("chr1", 25, 30)))
        self.assertTrue(obj.is_contact(CRange("chr1", 0, 10)))
        # Different chrom: always away, never overlap / adjacent / contact.
        other = CRange("chr2", 10, 20)
        self.assertFalse(obj.is_overlap(other))
        self.assertFalse(obj.is_adjacent(other))
        self.assertTrue(obj.is_away(other))
        self.assertFalse(obj.is_contact(other))
        # Even fully "overlapping" coordinates on different chromosomes.
        self.assertFalse(obj.is_overlap(CRange("chr2", 10, 20)))
        self.assertTrue(obj.is_away(CRange("chr2", 10, 20)))
        # Plain IRange operands fall back to the coordinate-only distance.
        self.assertTrue(obj.is_overlap(IRange(15, 25)))
        self.assertTrue(obj.is_away(IRange(25, 30)))

    def test_to_bed_string(self):
        """BED output: missing attributes (name/score/thick/color/blocks) render as None."""
        obj = CRange("chr1", 10, 20)
        # BED needs at least 3 columns.
        self.assertRaises(AssertionError, obj.to_bed_string, 1)
        self.assertEqual("chr1\t10\t20", obj.to_bed_string(ncol=3))
        # No name attribute: rendered as "None".
        self.assertEqual("chr1\t10\t20\tNone", obj.to_bed_string(ncol=4))
        # No score attribute: rendered as ".".
        self.assertEqual("chr1\t10\t20\tNone\t.", obj.to_bed_string(ncol=5))
        # No thick attribute: columns 7-8 fall back to start.
        self.assertEqual("chr1\t10\t20\tNone\t.\tNone\t10\t10", obj.to_bed_string(ncol=8))
        # No blocks attribute: rendered as a single block over the span.
        self.assertEqual(
            "chr1\t10\t20\tNone\t.\tNone\t10\t10\t255,0,0\t1\t10,\t0,",
            obj.to_bed_string(ncol=12))

    def test_to_gff_string(self):
        """GFF output: missing attributes (source/feature/score/frame/name) render as "."."""
        obj = CRange("chr1", 10, 20)
        # name is not set: printed as None in the attributes column.
        self.assertEqual(
            "chr1\t.\t.\t11\t20\t.\tNone\t.\tname=None",
            obj.to_gff_string())

    def test_to_tbl_string(self):
        """TBL output: plain attribute keys plus blocks/thick fallbacks."""
        obj = CRange("chr1", 10, 20)
        self.assertEqual("chr1\t10\t20", obj.to_tbl_string(header=["chrom", "start", "end"]))
        # No blocks attribute: the whole span is rendered as one block.
        self.assertEqual("10-20", obj.to_tbl_string(header=["blocks"]))
        # No thick attribute: rendered as an empty cell.
        self.assertEqual("", obj.to_tbl_string(header=["thick"]))
        # header is required.
        self.assertRaises(RuntimeError, obj.to_tbl_string, None)

    def test_format(self):
        """format() dispatches to the bed/gff/tbl writers."""
        obj = CRange("chr1", 10, 20)
        self.assertEqual("chr1\t10\t20", obj.format("bed", ncol=3))
        self.assertEqual(
            "chr1\t.\t.\t11\t20\t.\tNone\t.\tname=None",
            obj.format("gff"))
        self.assertEqual("chr1\t10\t20", obj.format("tbl", header=["chrom", "start", "end"]))
        # Format name is case-insensitive; unknown formats raise.
        self.assertEqual("chr1\t10\t20", obj.format("BED", ncol=3))
        self.assertRaises(RuntimeError, obj.format, "unknown")


if __name__ == "__main__":
    unittest.main()
