from pyBioInfo.Range import GRange


class Transcript(GRange):
    def __init__(self, obj=None, **kwargs):
        if obj is None:
            super(Transcript, self).__init__(**kwargs)
        else:
            super(Transcript, self).__init__(chrom=obj.chrom, start=obj.start, end=obj.end, name=obj.name,
                                             strand=obj.strand, blocks=obj.blocks, thick=obj.thick, **kwargs)

    @property
    def cds(self):
        return self.thick

    @property
    def exons(self):
        return self.exons

    @property
    def utr5(self):
        cds = self.cds
        if cds is None:
            return None
        if self.forward:
            return self.start, cds[0]
        else:
            return cds[1], self.end

    @property
    def utr3(self):
        cds = self.cds
        if cds is None:
            return None
        if self.forward:
            return cds[1], self.end
        else:
            return self.start, cds[0]

    def get_cds_length(self):
        cds = self.cds
        if cds is None:
            return None
        index1 = self.get_index(cds[0])
        index2 = self.get_index(cds[1] - 1)
        return index2 - index1 + 1

    def get_utr5_length(self):
        utr5 = self.utr5
        if utr5 is None:
            return None
        index1 = self.get_index(utr5[0])
        index2 = self.get_index(utr5[1] - 1)
        return index2 - index1 + 1

    def get_utr3_length(self):
        utr3 = self.utr3
        if utr3 is None:
            return None
        index1 = self.get_index(utr3[0])
        index2 = self.get_index(utr3[1] - 1)
        return index2 - index1 + 1

    def distance_base(self, position):
        i1 = 0
        i2 = None
        i3 = None
        i4 = len(self)
        cds = self.cds

        if cds is None:
            return None, None, None, None
        else:
            i2 = self.get_index(cds[0])
            i3 = self.get_index(cds[1] - 1)
            i5 = self.get_index(position)
            if self.forward:
                # i2 = i2
                i3 = i3 + 1
            else:
                i2, i3 = i3, i2 + 1
            return i5 - i1, i5 - i2, i5 - i3, i5 - i4

    def distance_ratio(self, position):
        i1 = 0
        i2 = None
        i3 = None
        i4 = len(self)
        cds = self.cds

        if cds is None:
            return None
        else:
            i2 = self.get_index(cds[0])
            i3 = self.get_index(cds[1] - 1)
            i5 = self.get_index(position)

            if self.forward:
                # i2 = i2
                i3 = i3 + 1
            else:
                i2, i3 = i3, i2 + 1
            ratio = None
            if i1 <= i5 < i2:
                ratio = (i5 - i1 + 0.5) / (i2 - i1)
            elif i2 <= i5 < i3:
                ratio = (i5 - i2 + 0.5) / (i3 - i2) + 1
            elif i3 <= i5 < i4:
                ratio = (i5 - i3 + 0.5) / (i4 - i3) + 2
            else:
                raise RuntimeError()
            return ratio
