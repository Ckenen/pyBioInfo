#!/usr/scripts/env python

import unittest
from pyBioInfo.Range import CRange


class TestChromosomeRange(unittest.TestCase):
    def test_chromsome_range(self):
        """
        Constructor
        """
        self.assertRaises(ValueError, CRange, chrom="chr1", start=10, end=20, strand="S")

        """
        Attributes
        """
        obj = CRange("chr1", 10, 20, "name1", "+")
        self.assertEqual("chr1", obj.chrom)
        self.assertEqual("name1", obj.name)
        self.assertEqual("+", obj.strand)
        self.assertTrue(obj.is_forward)
        self.assertFalse(obj.is_reverse)
        obj.chrom = "chr2"
        self.assertEqual("chr2", obj.chrom)
        obj.name = "name2"
        self.assertEqual("name2", obj.name)
        obj.strand = "-"
        self.assertEqual("-", obj.strand)
        self.assertFalse(obj.is_forward)
        self.assertTrue(obj.is_reverse)
        obj.strand = None
        self.assertIsNone(obj.strand)
        self.assertFalse(obj.is_forward)
        self.assertFalse(obj.is_reverse)
        
        """
        Index and position
        """
        obj = CRange("chr1", 10, 20, "name1", "+")
        array1 = list(range(10, 20))
        array2 = list(range(0, 10))
        array3 = list(range(-10, 0))
        for pos, idx in zip(array1, array2):
            self.assertEqual(idx, obj.get_index(pos))
        for idx, pos in zip(array2, array1):
            self.assertEqual(pos, obj.get_position(idx))
        for idx, pos in zip(array3, array1):
            self.assertEqual(pos, obj.get_position(idx))

        obj = CRange("chr1", 10, 20, "name1", "-")
        array1 = list(range(10, 20))
        array2 = list(range(0, 10))
        array3 = list(range(-10, 0))
        for pos, idx in zip(array1, array2[::-1]):
            self.assertEqual(idx, obj.get_index(pos))
        for idx, pos in zip(array2, array1[::-1]):
            self.assertEqual(pos, obj.get_position(idx))
        for idx, pos in zip(array3, array1[::-1]):
            self.assertEqual(pos, obj.get_position(idx))
        for pos, idx in zip(array1, array2):
            self.assertEqual(idx, obj.get_index(pos, strandness=False))
        for idx, pos in zip(array2, array1):
            self.assertEqual(pos, obj.get_position(idx, strandness=False))
        for idx, pos in zip(array3, array1):
            self.assertEqual(pos, obj.get_position(idx, strandness=False))


if __name__ == '__main__':
    unittest.main()
