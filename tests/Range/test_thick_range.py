# Testing module of pyBioInfo.Range.thick_range

import unittest

from pyBioInfo.Range import TRange


class TestTRange(unittest.TestCase):
    """Testing class of pyBioInfo.Range.thick_range.TRange."""

    def test_thick_attributes(self):
        """Thick interval storage and constructor validation."""
        blocks = [(10, 20), (30, 40)]
        # No thick by default.
        obj = TRange(chrom="chr1", blocks=blocks)
        self.assertIsNone(obj.thick)
        self.assertEqual((0, None, None, 20), obj.indexes())
        # Thick inside the blocks.
        obj = TRange(chrom="chr1", blocks=blocks, thick=(15, 35))
        self.assertEqual((15, 35), obj.thick)
        self.assertEqual((0, 5, 15, 20), obj.indexes())
        # Thick start/end out of the overall span.
        self.assertRaises(ValueError, TRange, chrom="chr1", blocks=blocks, thick=(5, 15))
        self.assertRaises(ValueError, TRange, chrom="chr1", blocks=blocks, thick=(35, 45))
        # Thick start/end not covered by a block.
        self.assertRaises(ValueError, TRange, chrom="chr1", blocks=blocks, thick=(20, 35))
        self.assertRaises(ValueError, TRange, chrom="chr1", blocks=blocks, thick=(15, 30))

    def test_thick_indexes_strandness(self):
        """Strand-aware thick indexes."""
        blocks = [(10, 20), (30, 40)]
        obj = TRange(chrom="chr1", blocks=blocks, thick=(12, 35), strand="-")
        self.assertEqual((0, 5, 18, 20), obj.indexes())
        self.assertEqual((0, 2, 15, 20), obj.indexes(strandness=False))

    def test_format_bed_thick(self):
        """BED thick columns filled when a thick interval is present."""
        obj = TRange(chrom="chr1", blocks=[(10, 20), (40, 50)], thick=(40, 50), name="name1", strand="+")
        self.assertEqual(
            "chr1\t10\t50\tname1\t.\t+\t40\t50\t0,0,255\t2\t10,10,\t0,30,",
            obj.format("bed"))
        # Without thick, columns 7-8 fall back to the overall start.
        obj2 = TRange(chrom="chr1", blocks=[(10, 20), (40, 50)], name="name1", strand="+")
        self.assertEqual(
            "chr1\t10\t50\tname1\t.\t+\t10\t10\t0,0,255\t2\t10,10,\t0,30,",
            obj2.format("bed"))

    def test_format_tbl(self):
        """TBL formatting with blocks and thick keys."""
        obj = TRange(chrom="chr1", blocks=[(10, 20), (40, 50)], thick=(40, 50))
        self.assertEqual(
            "chr1\t40-50\t10-20;40-50",
            obj.format("tbl", header=["chrom", "thick", "blocks"]))
        # No thick renders an empty cell.
        obj2 = TRange(chrom="chr1", blocks=[(10, 20)])
        self.assertEqual(
            "chr1\t\t10-20",
            obj2.format("tbl", header=["chrom", "thick", "blocks"]))
        # Plain interval keys.
        obj3 = TRange("chr1", 10, 20)
        self.assertEqual("chr1\t10\t20", obj3.format("tbl", header=["chrom", "start", "end"]))

    def test_string(self):
        """String representation includes the thick interval."""
        obj = TRange(chrom="chr1", blocks=[(10, 20), (40, 50)], thick=(40, 50), strand="+")
        self.assertEqual("TRange: 10-50 [strand: +, blocks: 2, thick: 40-50]", str(obj))
        # Without thick.
        obj2 = TRange(chrom="chr1", blocks=[(10, 20), (40, 50)])
        self.assertEqual("TRange: 10-50 [strand: None, blocks: 2]", str(obj2))


if __name__ == "__main__":
    unittest.main()
