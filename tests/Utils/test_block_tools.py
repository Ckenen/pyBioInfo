#!/usr/scripts/env python

import unittest
from pyBioInfo.Utils import BlockTools


class TestBlockTools(unittest.TestCase):
    def test_block_tools(self):
        """
        Length
        """
        blocks = [(10, 20), (30, 40), (50, 60)]
        self.assertEqual(BlockTools.get_length(blocks), 30)

        """
        Move
        """
        blocks1 = [(10, 20), (30, 40), (50, 60)]
        blocks2 = [(20, 30), (40, 50), (60, 70)]
        blocks3 = [(5, 15), (25, 35), (45, 55)]
        self.assertEqual(BlockTools.move(blocks1, 10), blocks2)
        self.assertEqual(BlockTools.move(blocks1, -5), blocks3)

        """
        Index and position
        """
        blocks = [(10, 20), (30, 40), (50, 60)]
        indexes1 = list(range(-30, 0))
        indexes2 = list(range(0, 30))
        positions = list(range(10, 20)) + list(range(30, 40)) + list(range(50, 60))
        for index1, index2, position in zip(indexes1, indexes2, positions):
            self.assertEqual(BlockTools.get_position(blocks, index1), position)
            self.assertEqual(BlockTools.get_position(blocks, index2), position)
            self.assertEqual(BlockTools.get_index(blocks, position), index2)
        self.assertRaises(ValueError, BlockTools.get_index, blocks, -1)
        self.assertRaises(ValueError, BlockTools.get_index, blocks, 9)
        self.assertRaises(ValueError, BlockTools.get_index, blocks, 21)
        self.assertRaises(ValueError, BlockTools.get_index, blocks, 29)
        self.assertRaises(ValueError, BlockTools.get_index, blocks, 41)
        self.assertRaises(ValueError, BlockTools.get_index, blocks, 49)
        self.assertRaises(ValueError, BlockTools.get_index, blocks, 61)
        self.assertRaises(ValueError, BlockTools.get_index, blocks, 100)
        self.assertRaises(ValueError, BlockTools.get_position, blocks, -31)
        self.assertRaises(ValueError, BlockTools.get_position, blocks, 30)

        """
        Gap
        """
        blocks1 = [(10, 20), (30, 40), (50, 60)]
        blocks2 = [(20, 30), (40, 50)]
        self.assertEqual(BlockTools.gaps(blocks1), blocks2)
        blocks1 = [(10, 20)]
        blocks2 = []
        self.assertEqual(BlockTools.gaps(blocks1), blocks2)
        blocks1 = []
        blocks2 = []
        self.assertEqual(BlockTools.gaps(blocks1), blocks2)

        """
        Suture
        """
        blocks1 = [(10, 15), (12, 20), (20, 30), (31, 40), (42, 50)]
        blocks2 = [(10, 30), (31, 40), (42, 50)]
        blocks3 = [(10, 40), (42, 50)]
        blocks4 = [(10, 50)]
        self.assertEqual(BlockTools.suture(blocks1, 0), blocks2)
        self.assertEqual(BlockTools.suture(blocks1, 1), blocks3)
        self.assertEqual(BlockTools.suture(blocks1, 2), blocks4)

        self.assertEqual([(5, 20), (30, 40), (50, 65)],
                         BlockTools.suture([(5, 10), (10, 20), (30, 40), (50, 60), (60, 65)]))

        """
        Fusion
        """
        blocks1 = [(10, 20), (30, 40), (50, 60)]
        blocks2 = [(15, 25), (40, 45), (65, 70)]
        blocks3 = [(10, 25), (30, 45), (50, 60), (65, 70)]
        self.assertEqual(BlockTools.fusion(blocks1, blocks2), blocks3)

        """
        Contain
        """
        blocks1 = [(10, 20), (40, 50), (60, 70)]
        blocks2 = [(10, 20), (40, 50), (60, 70)]
        blocks3 = [(15, 20), (40, 50), (60, 70)]
        blocks4 = [(15, 20), (42, 45), (60, 70)]
        blocks5 = [(10, 20), (40, 50), (60, 65)]
        blocks6 = [(8, 15), (45, 50), (60, 65)]
        blocks7 = [(10, 20), (45, 55), (60, 70)]
        self.assertTrue(BlockTools.is_contain(blocks1, blocks2))
        self.assertTrue(BlockTools.is_contain(blocks1, blocks3))
        self.assertTrue(BlockTools.is_contain(blocks1, blocks4))
        self.assertTrue(BlockTools.is_contain(blocks1, blocks5))
        self.assertFalse(BlockTools.is_contain(blocks1, blocks6))
        self.assertFalse(BlockTools.is_contain(blocks1, blocks7))

        """
        Coincide
        """
        blocks1 = [(10, 20), (40, 50), (60, 70)]
        blocks2 = [(10, 20), (40, 50), (60, 65)]
        blocks3 = [(15, 20), (40, 50), (60, 70)]
        blocks4 = [(5, 20), (40, 50), (60, 70)]
        blocks5 = [(10, 20), (40, 50), (60, 75)]
        blocks6 = [(10, 20), (40, 50), (62, 65)]
        blocks7 = [(10, 20), (40, 50), (58, 65)]
        blocks8 = [(10, 20), (40, 52), (60, 65)]
        blocks9 = [(10, 20), (40, 48), (62, 65)]
        blocks10 = [(10, 20), (40, 45), (45, 50), (60, 65)]
        self.assertTrue(BlockTools.is_coincide(blocks1, blocks2))
        self.assertTrue(BlockTools.is_coincide(blocks1, blocks3))
        self.assertFalse(BlockTools.is_coincide(blocks1, blocks4))
        self.assertFalse(BlockTools.is_coincide(blocks1, blocks5))
        self.assertFalse(BlockTools.is_coincide(blocks1, blocks6))
        self.assertFalse(BlockTools.is_coincide(blocks1, blocks7))
        self.assertFalse(BlockTools.is_coincide(blocks1, blocks8))
        self.assertFalse(BlockTools.is_coincide(blocks1, blocks9))
        self.assertFalse(BlockTools.is_coincide(blocks1, blocks10))
        blocks1 = [(10, 20)]
        blocks2 = [(12, 18)]
        blocks3 = [(8, 20)]
        blocks4 = [(10, 22)]
        self.assertTrue(BlockTools.is_coincide(blocks1, blocks2))
        self.assertFalse(BlockTools.is_coincide(blocks1, blocks3))
        self.assertFalse(BlockTools.is_coincide(blocks1, blocks4))

        """
        Clip
        """
        blocks1 = [(10, 20), (30, 40), (50, 60)]
        blocks2 = [(15, 20), (30, 40), (50, 55)]
        blocks3 = [(5, 20), (30, 40), (50, 65)]
        self.assertEqual(blocks2, BlockTools.clip(blocks1, 15, 55))
        self.assertEqual(blocks1, BlockTools.clip(blocks1, 5, 65))
        self.assertEqual(blocks3, BlockTools.clip(blocks1, 5, 65, extend=True, suture=True))

        blocks = [(10, 20), (30, 40), (50, 60)]
        template = [(0, 20), (30, 40), (50, 60), (65, 80)]
        expected = [(5, 20), (30, 40), (50, 60), (65, 70)]
        self.assertEqual(expected, BlockTools.clip(blocks, 5, 70, template=template, extend=True, suture=True))

        blocks = [(10, 20), (30, 40), (50, 60)]
        template = [(0, 20), (30, 40), (50, 60), (65, 80)]
        expected = [(5, 20), (30, 40), (50, 60), (65, 90)]
        self.assertEqual(expected,
                         BlockTools.clip(blocks, 5, 90, template=template, extend=True, supply=True, suture=True))

        """
        Extend
        """

        """
        Is sorted
        """
        blocks = [(10, 20), (30, 40), (50, 60)]
        self.assertTrue(BlockTools.is_sorted(blocks))
        blocks = [(10, 20), (50, 60), (30, 40)]
        self.assertFalse(BlockTools.is_sorted(blocks))

        """
        Check sorted
        """
        blocks = [(10, 20), (30, 40), (50, 60)]
        BlockTools.is_sorted(blocks)
        blocks = [(10, 20), (50, 60), (30, 40)]
        self.assertRaises(ValueError, BlockTools.check_sorted, blocks)

        """
        Is valid blocks
        """
        blocks = [(10, 20), (30, 40), (50, 60)]
        self.assertTrue(BlockTools.is_valid_blocks(blocks))
        blocks = [(10, 20), (15, 30)]
        self.assertFalse(BlockTools.is_valid_blocks(blocks))
        blocks = [(-5, 10), (20, 30)]
        self.assertFalse(BlockTools.is_valid_blocks(blocks))
        blocks = [(20, 10)]
        self.assertFalse(BlockTools.is_valid_blocks(blocks))

        """
        Check valid blocks
        """
        blocks = [(10, 20), (30, 40), (50, 60)]
        BlockTools.check_blocks(blocks)
        blocks = [(10, 20), (15, 30)]
        self.assertRaises(ValueError, BlockTools.check_blocks, blocks)
        blocks = [(-5, 10), (20, 30)]
        self.assertRaises(ValueError, BlockTools.check_blocks, blocks)
        blocks = [(20, 10)]
        self.assertRaises(ValueError, BlockTools.check_blocks, blocks)


if __name__ == '__main__':
    unittest.main()
