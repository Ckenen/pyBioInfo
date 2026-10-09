#!/usr/scripts/env python

import unittest
from pyBioInfo.Range import GRange


class TestGenomicRange(unittest.TestCase):
    def test_constructor(self):
        """Constructor validation."""
        blocks = [(10, 20), (40, 50), (60, 70)]
        # Inconsistent start/end with blocks.
        self.assertRaises(ValueError, GRange, chrom="chr1", start=10, end=60, blocks=blocks)
        self.assertRaises(ValueError, GRange, chrom="chr1", start=10, end=71, blocks=blocks)
        # A valid start of 0 must not skip the consistency check.
        self.assertRaises(ValueError, GRange, chrom="chr1", start=5, end=10, blocks=[(0, 10)])
        self.assertIsInstance(GRange(chrom="chr1", start=0, end=10, blocks=[(0, 10)]), GRange)
        # Invalid blocks.
        self.assertRaises(ValueError, GRange, chrom="chr1", blocks=[(20, 10)])
        self.assertRaises(ValueError, GRange, chrom="chr1", blocks=[(10, 20), (15, 25)])
        # Missing start/end without blocks.
        self.assertRaises(AssertionError, GRange, chrom="chr1")
        # Invalid strand (validated by the setter).
        self.assertRaises(ValueError, GRange, chrom="chr1", start=10, end=20, strand="S")
        self.assertRaises(ValueError, GRange, "chr1", -1, 20)  # start < 0
        self.assertRaises(ValueError, GRange, "chr1", 20, 10)  # start >= end
        self.assertRaises(TypeError, GRange, "chr1", "10", 20)  # non-integer start
        # Valid strands.
        for strand in (None, ".", "+", "-"):
            self.assertIsInstance(GRange("chr1", 10, 20, strand=strand), GRange)
        # Invalid strand via the setter after construction.
        obj = GRange("chr1", 10, 20)
        self.assertRaises(ValueError, setattr, obj, "strand", "S")

    def test_attributes(self):
        """Span, blocks, cached length and chromosome attributes."""
        blocks = [(10, 20), (40, 50), (60, 70)]
        obj = GRange(chrom="chr1", blocks=blocks)
        self.assertEqual(10, obj.start)
        self.assertEqual(70, obj.end)
        self.assertEqual(3, len(obj.blocks))
        self.assertEqual(30, len(obj))
        obj = GRange(chrom="chr1", start=10, end=20)
        self.assertEqual(10, len(obj))
        self.assertEqual(((10, 20),), obj.blocks)
        obj = GRange(chrom="chr1", start=10, end=70, blocks=blocks, name="name", strand="-")
        for block1, block2 in zip(blocks, obj.blocks):
            self.assertEqual(block1, block2)
        self.assertEqual(len(obj), 30)
        self.assertEqual(obj.chrom, "chr1")
        self.assertEqual(obj.name, "name")
        self.assertEqual(obj.strand, "-")
        self.assertTrue(obj.is_reverse)

    def test_chromosome_attributes(self):
        """Chromosome, name and strand accessors (merged from CRange)."""
        obj = GRange("chr1", 10, 20, "name1", "+")
        self.assertEqual(obj.chrom, "chr1")
        self.assertEqual(obj.name, "name1")
        self.assertEqual(obj.strand, "+")
        self.assertTrue(obj.is_forward)
        self.assertFalse(obj.is_reverse)
        obj.chrom = "chr2"
        self.assertEqual(obj.chrom, "chr2")
        obj.name = "name2"
        self.assertEqual(obj.name, "name2")
        obj.strand = "-"
        self.assertEqual(obj.strand, "-")
        self.assertFalse(obj.is_forward)
        self.assertTrue(obj.is_reverse)
        obj.strand = None
        self.assertIsNone(obj.strand)
        self.assertFalse(obj.is_forward)
        self.assertFalse(obj.is_reverse)

    def test_thick_attributes(self):
        """GRange no longer carries the thick attribute (moved to TRange)."""
        obj = GRange(chrom="chr1", blocks=[(10, 20), (30, 40)])
        self.assertFalse(hasattr(obj, "thick"))
        self.assertRaises(TypeError, GRange, chrom="chr1", blocks=[(10, 20), (30, 40)], thick=(15, 35))

    def test_format_tbl(self):
        """TBL formatting with plain interval keys."""
        obj = GRange("chr1", 10, 20)
        self.assertEqual("chr1\t10\t20", obj.format("tbl", header=["chrom", "start", "end"]))

    def test_sequence_protocol(self):
        """Iterable and indexable as a 2-element (start, end) sequence (merged from IRange)."""
        obj = GRange("chr1", 10, 20)
        self.assertEqual(10, obj.start)
        self.assertEqual(20, obj.end)
        self.assertEqual(10, obj[0])
        self.assertEqual(20, obj[1])
        start, end = obj
        self.assertEqual((10, 20), (start, end))
        self.assertEqual([10, 20], list(iter(obj)))
        self.assertRaises(RuntimeError, obj.__getitem__, 2)
        # Coordinates are read-only.
        self.assertRaises(AttributeError, setattr, obj, "start", 15)
        self.assertRaises(AttributeError, setattr, obj, "end", 25)
        self.assertRaises(RuntimeError, obj.__setitem__, 0, 15)
        self.assertRaises(RuntimeError, obj.__setitem__, 1, 25)

    def test_distance_relationships(self):
        """Distance relations between ranges (merged from IRange)."""
        obj = GRange("chr1", 10, 20)
        # (start, end): (is_overlap, is_adjacent, is_away, is_contact)
        matrix = [
            ((0, 5), (False, False, True, False)),    # strictly away
            ((0, 10), (False, True, False, True)),    # adjacent at the start
            ((0, 15), (True, False, False, True)),    # overlapping the start
            ((0, 25), (True, False, False, True)),    # covering
            ((5, 15), (True, False, False, True)),
            ((15, 25), (True, False, False, True)),   # overlapping the end
            ((20, 25), (False, True, False, True)),   # adjacent at the end
            ((25, 30), (False, False, True, False)),  # strictly away
        ]
        for (start, end), expected in matrix:
            other = GRange("chr1", start, end)
            overlap, adjacent, away, contact = expected
            self.assertEqual(overlap, obj.is_overlap(other), (start, end))
            self.assertEqual(adjacent, obj.is_adjacent(other), (start, end))
            self.assertEqual(away, obj.is_away(other), (start, end))
            self.assertEqual(contact, obj.is_contact(other), (start, end))
        # Sign convention of get_distance: max(starts) - min(ends), i.e.
        # negative = overlap length, 0 = adjacent, positive = gap length.
        self.assertEqual(-5, obj.get_distance(GRange("chr1", 15, 25)))
        self.assertEqual(0, obj.get_distance(GRange("chr1", 0, 10)))
        self.assertEqual(5, obj.get_distance(GRange("chr1", 25, 30)))

    def test_contain_single_block(self):
        """Coordinate containment for single-block ranges (merged from IRange)."""
        outer = GRange("chr1", 10, 30)
        self.assertTrue(outer.is_contain(GRange("chr1", 15, 25)))
        self.assertFalse(outer.is_contain(GRange("chr1", 5, 25)))

    def test_string(self):
        """String representation format."""
        obj = GRange(chrom="chr1", blocks=[(10, 20), (40, 50)], strand="+")
        self.assertEqual("GRange: 10-50 [strand: +, blocks: 2]", str(obj))

    def test_index_and_position(self):
        """Block-aware index/position conversion, gaps skipped."""
        blocks = [(10, 20), (40, 50), (60, 70)]
        obj = GRange(chrom="chr1", blocks=blocks, name="name", strand="-")
        array1 = list(range(10, 20)) + list(range(40, 50)) + list(range(60, 70))
        array2 = list(range(30))
        array3 = list(range(-30, 0))
        for pos, idx in zip(array1, array2):
            self.assertEqual(pos, obj.get_position(idx, strandness=False))
            self.assertEqual(idx, obj.get_index(pos, strandness=False))
        for pos, idx in zip(array1, array3):
            self.assertEqual(pos, obj.get_position(idx, strandness=False))

        for pos, idx in zip(array1[::-1], array2):
            self.assertEqual(pos, obj.get_position(idx, strandness=True))
            self.assertEqual(idx, obj.get_index(pos, strandness=True))
        for pos, idx in zip(array1[::-1], array3):
            self.assertEqual(pos, obj.get_position(idx, strandness=True))

        for position in (9, 21, 25, 71):  # in a gap or out of the span
            self.assertRaises(ValueError, obj.get_index, position)
        self.assertRaises(ValueError, obj.get_position, -31)
        self.assertRaises(ValueError, obj.get_position, 30)

    def test_index_and_position_forward(self):
        """Strand-aware index/position on the forward strand (single block)."""
        obj = GRange("chr1", 10, 20, "name1", "+")
        positions: list[int] = list(range(10, 20))
        indexes: list[int] = list(range(0, 10))
        negative: list[int] = list(range(-10, 0))
        for pos, idx in zip(positions, indexes):
            self.assertEqual(obj.get_index(pos), idx)
        for idx, pos in zip(indexes, positions):
            self.assertEqual(obj.get_position(idx), pos)
        for idx, pos in zip(negative, positions):
            self.assertEqual(obj.get_position(idx), pos)

    def test_index_and_position_reverse(self):
        """Strand-aware index/position on the reverse strand (single block)."""
        obj = GRange("chr1", 10, 20, "name1", "-")
        positions: list[int] = list(range(10, 20))
        indexes: list[int] = list(range(0, 10))
        negative: list[int] = list(range(-10, 0))
        # With strandness, the indexing runs from the end towards the start.
        for pos, idx in zip(positions, indexes[::-1]):
            self.assertEqual(obj.get_index(pos), idx)
        for idx, pos in zip(indexes, positions[::-1]):
            self.assertEqual(obj.get_position(idx), pos)
        for idx, pos in zip(negative, positions[::-1]):
            self.assertEqual(obj.get_position(idx), pos)
        # Without strandness, indexing behaves like the forward strand.
        for pos, idx in zip(positions, indexes):
            self.assertEqual(obj.get_index(pos, strandness=False), idx)
        for idx, pos in zip(indexes, positions):
            self.assertEqual(obj.get_position(idx, strandness=False), pos)
        for idx, pos in zip(negative, positions):
            self.assertEqual(obj.get_position(idx, strandness=False), pos)

    def test_comparison(self):
        """Ordering groups by chromosome first, then by coordinates."""
        obj1 = GRange("chr1", 10, 20)
        obj2 = GRange("chr1", 15, 25)
        obj3 = GRange("chr2", 5, 15)
        self.assertLess(obj1, obj2)
        self.assertLess(obj2, obj3)
        self.assertLess(obj1, obj3)
        self.assertGreater(obj3, obj2)
        self.assertEqual(obj1, GRange("chr1", 10, 20))
        self.assertNotEqual(obj1, obj2)
        # Same chrom and start: ordered by the end coordinate.
        self.assertLess(GRange("chr1", 10, 20), GRange("chr1", 10, 30))
        self.assertLessEqual(GRange("chr1", 10, 20), GRange("chr1", 10, 20))
        self.assertGreaterEqual(GRange("chr1", 10, 30), GRange("chr1", 10, 20))
        # Comparison with a non-range object.
        self.assertNotEqual(obj1, 123)

    def test_relationship(self):
        """Containment and coincidence against block and plain ranges."""
        obj1 = GRange(chrom="chr1", blocks=[(10, 20), (40, 50), (60, 70)])
        # Block-level containment: the middle block of obj1 is a gap for obj2.
        obj2 = GRange(chrom="chr1", blocks=[(10, 20), (60, 70)])
        self.assertTrue(obj1.is_contain(obj2))
        self.assertFalse(obj2.is_contain(obj1))
        # Coincide: block boundaries must align, allowing the first block
        # to start later, the last block to end earlier, and intermediate
        # block boundaries to differ by up to tolerance. The check is
        # directional and requires a consecutive run (no reference block
        # may be skipped).
        self.assertTrue(obj1.is_coincide(GRange(chrom="chr1", blocks=[(15, 20), (40, 50), (60, 65)])))
        self.assertFalse(obj1.is_coincide(obj2))  # obj2 skips the middle reference block
        self.assertFalse(obj2.is_coincide(obj1))
        self.assertFalse(obj1.is_coincide(GRange(chrom="chr1", blocks=[(10, 20), (40, 70)])))  # crosses the gap
        # Plain (non-block) ranges fall back to coordinate comparison.
        self.assertTrue(obj1.is_contain(GRange("chr1", 45, 50)))
        self.assertFalse(obj1.is_contain(GRange("chr1", 5, 15)))
        #self.assertTrue(obj1.is_coincide(GRange("chr1", 10, 70)))

    def test_format_bed(self):
        """BED formatting with different column counts."""
        obj = GRange(chrom="chr1", blocks=[(10, 20), (40, 50)], name="name1", strand="+")
        self.assertEqual(
            "chr1\t10\t50\tname1\t.\t+\t10\t10\t0,0,255\t2\t10,10,\t0,30,",
            obj.format("bed"))
        self.assertEqual(obj.format("BED"), obj.format("bed"))
        # 3-column BED keeps chrom/start/end only.
        self.assertEqual("chr1\t10\t50", obj.format("bed", ncol=3))
        # Reverse strand defaults to the red color.
        obj2 = GRange(chrom="chr1", blocks=[(10, 20)], name="n", strand="-")
        self.assertEqual("chr1\t10\t20\tn\t.\t-\t10\t10\t255,0,0\t1\t10,\t0,", obj2.format("bed"))
        # 6-column BED of a single-block range (merged from CRange tests).
        obj3 = GRange("chr1", 10, 20, "name1", "+")
        self.assertEqual(obj3.format("bed", ncol=6), "chr1\t10\t20\tname1\t.\t+")
        self.assertEqual(GRange("chr1", 10, 20).format("bed", ncol=6), "chr1\t10\t20\tNone\t.\tNone")

    def test_format_gff(self):
        """GFF formatting."""
        obj = GRange(chrom="chr1", blocks=[(10, 20), (40, 50)], name="name1", strand="+")
        # 1-based start, missing annotation columns rendered as ".".
        self.assertEqual(
            "chr1\t.\t.\t11\t50\t.\t+\t.\tname=name1",
            obj.format("gff"))

    def test_format_others(self):
        """Unsupported formats and missing headers."""
        obj = GRange(chrom="chr1", blocks=[(10, 20)])
        self.assertRaises(RuntimeError, obj.format, "sam")
        self.assertRaises(RuntimeError, obj.format, "tbl")  # header required


if __name__ == '__main__':
    _ = unittest.main()
