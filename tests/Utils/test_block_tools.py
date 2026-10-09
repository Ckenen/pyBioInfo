#!/usr/scripts/env python

import unittest
from pyBioInfo.Utils import BlockTools


class TestBlockTools(unittest.TestCase):
    def test_length(self):
        """Summation length of all blocks."""
        blocks = [(10, 20), (30, 40), (50, 60)]
        self.assertEqual(BlockTools.get_length(blocks), 30)
        self.assertEqual(BlockTools.get_length([]), 0)

    def test_move(self):
        """Moving blocks by a signed step."""
        blocks1 = [(10, 20), (30, 40), (50, 60)]
        blocks2 = [(20, 30), (40, 50), (60, 70)]
        blocks3 = [(5, 15), (25, 35), (45, 55)]
        self.assertEqual(BlockTools.move(blocks1, 10), blocks2)
        self.assertEqual(BlockTools.move(blocks1, -5), blocks3)
        self.assertEqual(BlockTools.move(blocks1, 0), blocks1)
        self.assertRaises(RuntimeError, BlockTools.move, blocks1, -11)

    def test_index_and_position(self):
        """Block-aware index/position conversion, gaps skipped."""
        blocks = [(10, 20), (30, 40), (50, 60)]
        indexes1 = list(range(-30, 0))
        indexes2 = list(range(0, 30))
        positions = list(range(10, 20)) + list(range(30, 40)) + list(range(50, 60))
        for index1, index2, position in zip(indexes1, indexes2, positions):
            with self.subTest(position=position):
                self.assertEqual(BlockTools.get_position(blocks, index1), position)
                self.assertEqual(BlockTools.get_position(blocks, index2), position)
                self.assertEqual(BlockTools.get_index(blocks, position), index2)
        for position in (-1, 9, 21, 29, 41, 49, 61, 100):
            self.assertRaises(ValueError, BlockTools.get_index, blocks, position)
        self.assertRaises(ValueError, BlockTools.get_position, blocks, -31)
        self.assertRaises(ValueError, BlockTools.get_position, blocks, 30)
        # Reverse strand indexes run from the last block end backwards.
        self.assertEqual(BlockTools.get_index(blocks, 59, reverse=True), 0)
        self.assertEqual(BlockTools.get_position(blocks, 0, reverse=True), 59)

    def test_is_include(self):
        """Position membership test."""
        blocks = [(10, 20), (30, 40), (50, 60)]
        for position in (10, 15, 19, 30, 39, 50, 60 - 1):
            self.assertTrue(BlockTools.is_include(blocks, position))
        for position in (-1, 9, 20, 29, 40, 49, 60, 100):
            self.assertFalse(BlockTools.is_include(blocks, position))

    def test_gaps(self):
        """Gaps between adjacent blocks."""
        blocks1 = [(10, 20), (30, 40), (50, 60)]
        blocks2 = [(20, 30), (40, 50)]
        self.assertEqual(BlockTools.gaps(blocks1), blocks2)
        self.assertEqual(BlockTools.gaps([(10, 20)]), [])
        self.assertEqual(BlockTools.gaps([]), [])

    def test_suture(self):
        """Merging blocks whose gap is small enough."""
        blocks1 = [(10, 15), (12, 20), (20, 30), (31, 40), (42, 50)]
        self.assertEqual(BlockTools.suture(blocks1, 0), [(10, 30), (31, 40), (42, 50)])
        self.assertEqual(BlockTools.suture(blocks1, 1), [(10, 40), (42, 50)])
        self.assertEqual(BlockTools.suture(blocks1, 2), [(10, 50)])
        self.assertEqual(BlockTools.suture([]), [])
        self.assertEqual([(5, 20), (30, 40), (50, 65)],
                         BlockTools.suture([(5, 10), (10, 20), (30, 40), (50, 60), (60, 65)]))

    def test_fusion(self):
        """Merging two block lists into one."""
        blocks1 = [(10, 20), (30, 40), (50, 60)]
        blocks2 = [(15, 25), (40, 45), (65, 70)]
        blocks3 = [(10, 25), (30, 45), (50, 60), (65, 70)]
        self.assertEqual(BlockTools.fusion(blocks1, blocks2), blocks3)

    def test_contain(self):
        """Full coverage test."""
        blocks1 = [(10, 20), (40, 50), (60, 70)]
        contained = [
            [(10, 20), (40, 50), (60, 70)],
            [(15, 20), (40, 50), (60, 70)],
            [(15, 20), (42, 45), (60, 70)],
            [(10, 20), (40, 50), (60, 65)],
        ]
        for blocks2 in contained:
            self.assertTrue(BlockTools.is_contain(blocks1, blocks2))
        not_contained = [
            [(8, 15), (45, 50), (60, 65)],  # left boundary exceeds
            [(10, 20), (45, 55), (60, 70)],  # right boundary exceeds
        ]
        for blocks2 in not_contained:
            self.assertFalse(BlockTools.is_contain(blocks1, blocks2))

    def test_coincide(self):
        """Consecutive-run coincidence with boundary tolerance."""
        blocks1 = [(10, 20), (40, 50), (60, 70)]
        coincide = [
            [(10, 20), (40, 50), (60, 70)],   # identical
            [(15, 20), (40, 50), (60, 65)],   # first start / last end inset
            [(10, 20), (40, 50)],             # subset, consecutive run
            [(40, 50), (60, 70)],             # starts at the 2nd reference block
        ]
        for blocks2 in coincide:
            self.assertTrue(BlockTools.is_coincide(blocks1, blocks2))
        not_coincide = [
            [(10, 20), (40, 52), (60, 65)],   # middle block extends into the gap
            [(10, 20), (40, 48), (62, 65)],   # double-sided mismatch
            [(10, 20), (40, 45), (45, 50), (60, 65)],  # extra block in the gap
            [(10, 20), (50, 60)],             # skips a reference block
            [(10, 20), (40, 50), (60, 70), (80, 90)],  # more blocks than reference
        ]
        for blocks2 in not_coincide:
            self.assertFalse(BlockTools.is_coincide(blocks1, blocks2))
        # Single-block cases.
        self.assertTrue(BlockTools.is_coincide([(10, 20)], [(12, 18)]))
        self.assertFalse(BlockTools.is_coincide([(10, 20)], [(20, 30)]))
        self.assertFalse(BlockTools.is_coincide([(10, 20)], [(5, 8)]))
        # Tolerance allows intermediate boundary deviation.
        self.assertTrue(BlockTools.is_coincide(blocks1, [(10, 20), (42, 50), (60, 70)], tolerance=2))
        self.assertFalse(BlockTools.is_coincide(blocks1, [(10, 20), (42, 50), (60, 70)], tolerance=1))

    def test_boundary_bias(self):
        """Boundary bias of blocks_que relative to blocks_ref (signed)."""
        ref = [(10, 20), (30, 40), (50, 60)]
        # Start outside on the left -> negative coordinate difference.
        self.assertEqual(BlockTools.get_boundary_bias(ref, [(5, 20), (30, 40), (50, 60)]), (-5, 0))
        # Start inside a block -> positive block-aware index.
        self.assertEqual(BlockTools.get_boundary_bias(ref, [(15, 20), (30, 40), (50, 60)]), (5, 0))
        self.assertEqual(BlockTools.get_boundary_bias(ref, [(35, 40), (50, 60)]), (15, 0))
        # End outside on the right -> positive coordinate difference.
        self.assertEqual(BlockTools.get_boundary_bias(ref, [(10, 20), (30, 40), (50, 65)]), (0, 5))
        # End inside a block -> negative block-aware index.
        self.assertEqual(BlockTools.get_boundary_bias(ref, [(10, 20), (30, 40), (50, 55)]), (0, -6))
        # Identical boundaries.
        self.assertEqual(BlockTools.get_boundary_bias(ref, ref), (0, 0))
        # A query boundary falling inside a gap of the reference raises.
        self.assertRaises(ValueError, BlockTools.get_boundary_bias, ref, [(25, 40), (50, 60)])
        self.assertRaises(ValueError, BlockTools.get_boundary_bias, ref, [(10, 20), (25, 28)])
        # Invalid block lists are rejected.
        self.assertRaises(ValueError, BlockTools.get_boundary_bias, ref, [(10, 20), (15, 30)])

    def test_is_matched_gap(self):
        """Gaps of blocks_que match an ordered subset of blocks_ref gaps."""
        ref = [(10, 20), (30, 40), (50, 60)]
        # Identical gaps.
        self.assertTrue(BlockTools.is_matched_gap(ref, [(10, 20), (30, 40), (50, 60)]))
        # Fewer blocks, gaps are an ordered subset.
        self.assertTrue(BlockTools.is_matched_gap(ref, [(10, 20), (30, 40)]))
        self.assertTrue(BlockTools.is_matched_gap(ref, [(30, 40), (50, 60)]))
        # A gap spanning two reference gaps is not a subset.
        self.assertFalse(BlockTools.is_matched_gap(ref, [(10, 20), (50, 60)]))
        # The first query gap may match a later reference gap (not the first).
        ref4 = [(10, 20), (30, 40), (50, 60), (70, 80)]
        self.assertTrue(BlockTools.is_matched_gap(ref4, [(30, 40), (50, 60), (70, 80)]))
        self.assertTrue(BlockTools.is_matched_gap(ref4, [(50, 60), (70, 80)]))
        # Once the first gap matches, the rest must match consecutively.
        ref5 = [(10, 20), (30, 40), (50, 60), (70, 80), (90, 100)]
        self.assertFalse(BlockTools.is_matched_gap(ref5, [(30, 40), (50, 60), (90, 100)]))
        # Single block -> no gaps -> always matched.
        self.assertTrue(BlockTools.is_matched_gap(ref, [(10, 20)]))
        # More blocks than the reference -> False.
        self.assertFalse(BlockTools.is_matched_gap(ref, [(10, 20), (30, 40), (50, 60), (70, 80)]))
        # A gap that is not a subset of the reference gaps -> False.
        self.assertFalse(BlockTools.is_matched_gap(ref, [(10, 20), (25, 40), (50, 60)]))
        # Tolerance allows boundary deviation.
        self.assertTrue(BlockTools.is_matched_gap(ref, [(10, 20), (33, 40), (50, 60)], tolerance=3))
        self.assertFalse(BlockTools.is_matched_gap(ref, [(10, 20), (33, 40), (50, 60)], tolerance=2))
        # Invalid blocks are rejected.
        self.assertRaises(ValueError, BlockTools.is_matched_gap, ref, [(10, 20), (15, 30)])

    def test_clip(self):
        """Clipping blocks to a window, with optional extension/supply."""
        blocks1 = [(10, 20), (30, 40), (50, 60)]
        blocks2 = [(15, 20), (30, 40), (50, 55)]
        blocks3 = [(5, 20), (30, 40), (50, 65)]
        self.assertEqual(blocks2, BlockTools.clip(blocks1, 15, 55))
        self.assertEqual(blocks1, BlockTools.clip(blocks1, 5, 65))
        self.assertEqual(blocks3, BlockTools.clip(blocks1, 5, 65, extend=True, suture=True))

        template = [(0, 20), (30, 40), (50, 60), (65, 80)]
        expected = [(5, 20), (30, 40), (50, 60), (65, 70)]
        self.assertEqual(expected, BlockTools.clip(blocks1, 5, 70, template=template, extend=True, suture=True))

        expected = [(5, 20), (30, 40), (50, 60), (65, 90)]
        self.assertEqual(expected,
                         BlockTools.clip(blocks1, 5, 90, template=template, extend=True, supply=True, suture=True))

        # Empty window or negative coordinates are invalid.
        self.assertEqual(BlockTools.clip(blocks1, 20, 30), [])
        self.assertRaises(ValueError, BlockTools.clip, blocks1, -1, 30)
        self.assertRaises(ValueError, BlockTools.clip, blocks1, 30, 30)

    def test_sorting(self):
        """Sorted check and block sorting."""
        sorted_blocks = [(10, 20), (30, 40), (50, 60)]
        self.assertTrue(BlockTools.is_sorted(sorted_blocks))
        self.assertFalse(BlockTools.is_sorted([(10, 20), (50, 60), (30, 40)]))
        self.assertIsNone(BlockTools.check_sorted(sorted_blocks))
        self.assertRaises(ValueError, BlockTools.check_sorted, [(10, 20), (50, 60), (30, 40)])
        self.assertEqual(BlockTools.sorted([(30, 40), (10, 20), (50, 60)]), sorted_blocks)
        # Blocks with the same start are ordered by end.
        self.assertEqual(BlockTools.sorted([(10, 25), (10, 20), (10, 22)]),
                         [(10, 20), (10, 22), (10, 25)])

    def test_valid_blocks(self):
        """Block validity check and its boolean counterpart."""
        valid = [(10, 20), (30, 40), (50, 60)]
        invalid = [
            [(10, 20), (15, 30)],  # overlapping
            [(-5, 10), (20, 30)],  # negative start
            [(20, 10)],            # start >= end
        ]
        self.assertTrue(BlockTools.is_valid_blocks(valid))
        self.assertFalse(BlockTools.is_valid_blocks([]))
        self.assertIsNone(BlockTools.check_blocks(valid))
        for blocks in invalid:
            self.assertFalse(BlockTools.is_valid_blocks(blocks))
            self.assertRaises(ValueError, BlockTools.check_blocks, blocks)
        self.assertRaises(RuntimeError, BlockTools.check_blocks, [])  # empty


if __name__ == '__main__':
    _ = unittest.main()
