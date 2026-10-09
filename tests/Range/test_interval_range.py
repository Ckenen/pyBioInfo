# Testing module of pyBioInfo.Range.interval_range

import unittest

from pyBioInfo.Range import IRange


class TestIRange(unittest.TestCase):
    """Testing class of pyBioInfo.Range.interval_range.IRange."""

    def test_constructor(self):
        """Construction with read-only start and end attributes."""
        obj = IRange(5, 15)
        self.assertEqual(5, obj.start)
        self.assertEqual(15, obj.end)
        # start must be >= 0 and start < end.
        self.assertRaises(ValueError, IRange, -1, 10)
        self.assertRaises(ValueError, IRange, 20, 10)
        self.assertRaises(ValueError, IRange, 10, 10)
        # Non-integer coordinates are rejected by the comparison.
        self.assertRaises(TypeError, IRange, "10", 20)

    def test_read_only_attributes(self):
        """The start and end attributes can not be overwritten."""
        obj = IRange(5, 15)
        self.assertRaises(AttributeError, setattr, obj, "start", 8)
        self.assertRaises(AttributeError, setattr, obj, "end", 18)

    def test_distance_relationships(self):
        """Distance relations and the get_distance sign convention."""
        obj = IRange(10, 20)
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
            other = IRange(start, end)
            overlap, adjacent, away, contact = expected
            self.assertEqual(overlap, obj.is_overlap(other), (start, end))
            self.assertEqual(adjacent, obj.is_adjacent(other), (start, end))
            self.assertEqual(away, obj.is_away(other), (start, end))
            self.assertEqual(contact, obj.is_contact(other), (start, end))
        # Sign convention of get_distance: max(starts) - min(ends), i.e.
        # negative = overlap length, 0 = adjacent, positive = gap length.
        self.assertEqual(-5, obj.get_distance(IRange(15, 25)))
        self.assertEqual(0, obj.get_distance(IRange(0, 10)))
        self.assertEqual(5, obj.get_distance(IRange(25, 30)))

    def test_index_and_position(self):
        """Index/position conversion between coordinates and offsets."""
        obj = IRange(10, 20)
        # Coordinate -> 0-based index.
        self.assertEqual(0, obj.get_index(10))
        self.assertEqual(2, obj.get_index(12))
        self.assertEqual(9, obj.get_index(19))
        self.assertRaises(ValueError, obj.get_index, 9)
        self.assertRaises(ValueError, obj.get_index, 20)
        # Index -> coordinate, negative indexes count from the end.
        self.assertEqual(10, obj.get_position(0))
        self.assertEqual(12, obj.get_position(2))
        self.assertEqual(19, obj.get_position(9))
        self.assertEqual(19, obj.get_position(-1))
        self.assertEqual(10, obj.get_position(-10))
        self.assertRaises(ValueError, obj.get_position, 10)
        self.assertRaises(ValueError, obj.get_position, -11)

    def test_sequence_protocol(self):
        """len(), str(), iteration and indexing as a 2-element sequence."""
        obj = IRange(10, 20)
        self.assertEqual(10, len(obj))
        self.assertEqual("IRange: 10-20", str(obj))
        start, end = obj
        self.assertEqual((10, 20), (start, end))
        self.assertEqual([10, 20], list(iter(obj)))
        self.assertEqual(10, obj[0])
        self.assertEqual(20, obj[1])
        self.assertRaises(RuntimeError, obj.__getitem__, 2)
        self.assertRaises(RuntimeError, obj.__setitem__, 0, 15)
        self.assertRaises(RuntimeError, obj.__setitem__, 1, 25)

    def test_comparison(self):
        """Ordering by start first, then by end (coordinates only)."""
        obj1 = IRange(10, 20)
        obj2 = IRange(15, 25)
        self.assertLess(obj1, obj2)
        self.assertGreater(obj2, obj1)
        self.assertLessEqual(obj1, IRange(10, 20))
        self.assertGreaterEqual(obj2, IRange(15, 25))
        # Same start: ordered by the end coordinate.
        self.assertLess(IRange(10, 20), IRange(10, 30))
        # Equality of the sort key.
        self.assertEqual(obj1, IRange(10, 20))
        self.assertNotEqual(obj1, obj2)
        self.assertNotEqual(obj1, 123)
        # Usable with sorted / min.
        items = [IRange(20, 30), IRange(5, 15), IRange(10, 25), IRange(5, 10)]
        self.assertEqual([(5, 10), (5, 15), (10, 25), (20, 30)], [tuple(x) for x in sorted(items)])
        self.assertEqual((5, 10), tuple(min(items)))


if __name__ == "__main__":
    unittest.main()
