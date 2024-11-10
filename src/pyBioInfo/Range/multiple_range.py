# Class: MRange (multiple genomic range()
from .interval_range import IRange
from .genomic_range import GRange
from pyBioInfo.Utils import BlockTools


class MRange(GRange):
    def __init__(self, chrom=None, start=None, end=None, name=None, strand=None, 
                 blocks=None, blocks_array=None):
        blocks = None
        if blocks_array is None:
            if start is None:
                raise ValueError()
            if end is None:
                raise ValueError()
            blocks = tuple([tuple([start, end])])
            blocks_array = tuple([tuple([tuple([start, end])])])
        else:
            array1 = []
            array2 = []
            for blocks1 in blocks_array:
                BlockTools.check_blocks(blocks1)
                for bstart, bend in blocks1:
                    array2.append([bstart, bend])
                blocks1 = tuple([(bstart, bend) for bstart, bend in blocks1])
                array1.append(blocks1)
            array2 = BlockTools.sorted(array2)
            blocks = tuple([(bstart, bend) for bstart, bend in BlockTools.suture(array2)])
        assert blocks
        super(MRange, self).__init__(chrom=chrom, start=start, end=end, 
                                     name=name, strand=strand, blocks=blocks)
        self._blocks_array = blocks_array

    def coincide_from(self, other):
        for blocks in self._blocks_array:
            if not BlockTools.is_coincide(other.blocks, blocks):
                return False
        return True
