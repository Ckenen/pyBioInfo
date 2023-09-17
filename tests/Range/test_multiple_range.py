#!/usr/scripts/env python
import unittest
from pyBioInfo.Range import MRange


class TestMultipleRange(unittest.TestCase):
    def test_multiple_range(self):
        """
        Constructor
        """
        blocks1 = [(10, 20), (30, 40), (50, 60)]
        blocks2 = [(15, 25), (35, 38), (45, 65)]
        blocks3 = [(15, 20), (30, 40), (80, 90)]
        array = [blocks1, blocks2, blocks3]
        obj = MRange(chrom="chr1", name="name", strand="+", blocks_array=array)
        blocks = [(10, 25), (30, 40), (45, 65), (80, 90)]
        for block1, block2 in zip(blocks, obj.blocks):
            self.assertEqual(block1[0], block2[0])
            self.assertEqual(block1[1], block2[1])
        self.assertEqual(55, len(obj))

        """
        Coincide
        """
        blocks = [(10, 20), (30, 40), (50, 60)]
        array = [blocks]
        obj1 = MRange(chrom="chr1", blocks_array=array)

        blocks1 = [(15, 20), (30, 35)]
        blocks2 = [(35, 40), (50, 55)]
        array = [blocks1, blocks2]
        obj2 = MRange(chrom="chr1", blocks_array=array)
        self.assertTrue(obj2.coincide_from(obj1))

        blocks1 = [(15, 20), (30, 40)]
        blocks2 = [(35, 40), (50, 55)]
        array = [blocks1, blocks2]
        obj2 = MRange(chrom="chr1", blocks_array=array)
        self.assertTrue(obj2.coincide_from(obj1))

        blocks1 = [(15, 20), (32, 40)]
        blocks2 = [(35, 40), (50, 55)]
        array = [blocks1, blocks2]
        obj2 = MRange(chrom="chr1", blocks_array=array)
        self.assertFalse(obj2.coincide_from(obj1))

        blocks1 = [(15, 20), (30, 35)]
        blocks2 = [(35, 40), (50, 55)]
        blocks3 = [(55, 58)]
        array = [blocks1, blocks2, blocks3]
        obj2 = MRange(chrom="chr1", blocks_array=array)
        self.assertTrue(obj2.coincide_from(obj1))

        blocks1 = [(15, 20), (25, 35)]
        blocks2 = [(35, 40), (50, 55)]
        array = [blocks1, blocks2]
        obj2 = MRange(chrom="chr1", blocks_array=array)
        self.assertFalse(obj2.coincide_from(obj1))


if __name__ == '__main__':
    unittest.main()
