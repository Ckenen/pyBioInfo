#!/usr/scripts/env python

import unittest
from pyBioInfo.Range import IRange


class TestIntervalRange(unittest.TestCase):
    def test_interval_range(self):
        """
        Constructor
        """
        self.assertRaises(ValueError, IRange, -1, 10)  # start < 0
        self.assertRaises(ValueError, IRange, 10, 9)  # start < end
        self.assertRaises(ValueError, IRange, 10, 10)  # start == end
        self.assertRaises(TypeError, IRange, "10", 20)

        """
        Attributes
        """
        obj = IRange(10, 20)
        self.assertEqual(len(obj), 10)
        # Getter
        self.assertEqual(obj.start, 10)  # start == 10
        self.assertEqual(obj.end, 20)  # end == 20
        self.assertEqual(len(obj), 20 - 10)  # len()
        self.assertEqual(obj[0], 10)
        self.assertEqual(obj[1], 20)
        self.assertEqual(tuple([value for value in obj]), (10, 20))
        # Setter
        self.assertRaises(AttributeError, obj.__setattr__, "start", 5)
        self.assertRaises(AttributeError, obj.__setattr__, "end", 25)
        self.assertRaises(AttributeError, obj.__setattr__, "name", "value")

        """
        Index and position
        """
        obj = IRange(10, 20)
        indexes1 = list(range(10))
        indexes2 = list(range(-10, 0))
        positions = list(range(10, 20))
        for index1, index2, position in zip(indexes1, indexes2, positions):
            self.assertEqual(obj.get_index(position), index1)
            self.assertEqual(obj.get_position(index1), position)
            self.assertEqual(obj.get_position(index2), position)
        self.assertRaises(ValueError, obj.get_index, 9)
        self.assertRaises(ValueError, obj.get_index, 21)
        self.assertRaises(ValueError, obj.get_position, -11)
        self.assertRaises(ValueError, obj.get_position, 10)

        """
        Comparison
        """
        obj1 = IRange(10, 20)
        obj2 = IRange(10, 25)
        obj3 = IRange(5, 25)
        obj4 = IRange(15, 20)
        self.assertGreater(obj1, obj3)
        self.assertGreater(obj2, obj1)
        self.assertGreater(obj2, obj3)
        self.assertGreater(obj4, obj1)
        self.assertGreater(obj4, obj2)
        self.assertGreater(obj4, obj3)
        self.assertLess(obj3, obj1)
        self.assertLess(obj3, obj2)
        self.assertLess(obj3, obj4)

        """
        Distance
        """
        obj1 = IRange(10, 20)
        obj2 = IRange(15, 25)
        obj3 = IRange(20, 30)
        obj4 = IRange(25, 30)
        # overlap
        self.assertTrue(obj1.is_overlap(obj2))
        self.assertFalse(obj1.is_overlap(obj3))
        self.assertFalse(obj1.is_overlap(obj4))
        self.assertTrue(obj2.is_overlap(obj1))
        self.assertTrue(obj2.is_overlap(obj3))
        self.assertFalse(obj2.is_overlap(obj4))
        self.assertFalse(obj3.is_overlap(obj1))
        self.assertTrue(obj3.is_overlap(obj2))
        self.assertTrue(obj3.is_overlap(obj4))
        self.assertFalse(obj4.is_overlap(obj1))
        self.assertFalse(obj4.is_overlap(obj2))
        self.assertTrue(obj4.is_overlap(obj3))
        # away
        self.assertFalse(obj1.is_away(obj2))
        self.assertFalse(obj1.is_away(obj3))
        self.assertTrue(obj1.is_away(obj4))
        self.assertFalse(obj2.is_away(obj1))
        self.assertFalse(obj2.is_away(obj3))
        self.assertFalse(obj2.is_away(obj4))
        self.assertFalse(obj3.is_away(obj1))
        self.assertFalse(obj3.is_away(obj2))
        self.assertFalse(obj3.is_away(obj4))
        self.assertTrue(obj4.is_away(obj1))
        self.assertFalse(obj4.is_away(obj2))
        self.assertFalse(obj4.is_away(obj3))
        # dajacent
        self.assertFalse(obj1.is_adjacent(obj2))
        self.assertTrue(obj1.is_adjacent(obj3))
        self.assertFalse(obj1.is_adjacent(obj4))
        self.assertFalse(obj2.is_adjacent(obj1))
        self.assertFalse(obj2.is_adjacent(obj3))
        self.assertTrue(obj2.is_adjacent(obj4))
        self.assertTrue(obj3.is_adjacent(obj1))
        self.assertFalse(obj3.is_adjacent(obj2))
        self.assertFalse(obj3.is_adjacent(obj4))
        self.assertFalse(obj4.is_adjacent(obj1))
        self.assertTrue(obj4.is_adjacent(obj2))
        self.assertFalse(obj4.is_adjacent(obj3))


if __name__ == '__main__':
    unittest.main()
