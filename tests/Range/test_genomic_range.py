import unittest
from pyBioInfo.Range import GRange


class TestGenomicRange(unittest.TestCase):
    def test_genomic_range(self):
        """
        Constructor
        """
        blocks = [(10, 20), (40, 50), (60, 70)]
        self.assertRaises(ValueError, GRange, chrom="chr1", start=10, end=60, blocks=blocks)

        """
        Attributes
        """
        blocks = [(10, 20), (40, 50), (60, 70)]
        obj = GRange(chrom="chr1", blocks=blocks)
        self.assertEqual(10, obj.start)
        self.assertEqual(70, obj.end)
        self.assertEqual(3, len(obj.blocks))
        self.assertEqual(30, len(obj))
        obj = GRange(chrom="chr1", start=10, end=20)
        self.assertEqual(10, len(obj))
        obj = GRange(chrom="chr1", start=10, end=70, blocks=blocks, name="name", strand="-")
        for block1, block2 in zip(blocks, obj.blocks):
            self.assertEqual(block1, block2)
        self.assertEqual(len(obj), 30)
        self.assertEqual(obj.chrom, "chr1")
        self.assertEqual(obj.name, "name")

        """
        Index and position
        """
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

        self.assertRaises(ValueError, obj.get_index, 9)
        self.assertRaises(ValueError, obj.get_index, 21)
        self.assertRaises(ValueError, obj.get_index, 25)
        self.assertRaises(ValueError, obj.get_index, 71)
        self.assertRaises(ValueError, obj.get_position, -31)
        self.assertRaises(ValueError, obj.get_position, 30)


if __name__ == '__main__':
    unittest.main()
