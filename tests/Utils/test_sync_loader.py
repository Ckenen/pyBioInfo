#!/usr/scripts/env python

import unittest
from pyBioInfo.Utils import SyncLoader


class TestSyncLoader(unittest.TestCase):
    def test_sync_loader(self):
        items1 = [1, 2, 3, 4, 5, 6, 10]
        items2 = [2, 3, 4, 5, 5, 5]
        items3 = [0, 4, 4, 9, 9]
        items4 = []
        array = [items1, items2, items3, items4]
            
        loader = SyncLoader(array)
        for i, (index, item) in enumerate(loader):
            if i > 100:
                break
            # print(i, index, item, sep="\t")


if __name__ == '__main__':
    unittest.main()


