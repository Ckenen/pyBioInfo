from pyBioInfo.Range import GRange
from pyBioInfo.Utils import BlockTools


class Fragment(GRange):
    def __init__(self, mates=None):
        chrom = mates[0].chrom
        name = mates[0].name
        strand = mates[0].strand

        assert len(mates) == 2
        blocks = BlockTools.fusion(mates[0].blocks, mates[1].blocks)

        super(Fragment, self).__init__(chrom=chrom, name=name, strand=strand, blocks=blocks)

        self._mates = tuple(mates)

    @property
    def mates(self):
        return self._mates

    # def blocks(self, *args, **kwargs):
    #     if self.is_single_end:
    #         for block in self.mate1.blocks():
    #             yield block
    #     elif self.is_paired_end:
    #         blocks1 = list(self.mate1.blocks())
    #         blocks2 = list(self.mate2.blocks())
    #         for block in BlockTools.fusion(blocks1, blocks2):
    #             yield block
    #     else:
    #         raise RuntimeError()

    @property
    def is_single_end(self):
        return len(self.mates) == 1

    @property
    def is_paired_end(self):
        return len(self.mates) == 2

    @property
    def mate1(self):
        return self.mates[0]

    @property
    def mate2(self):
        return self.mates[1]
