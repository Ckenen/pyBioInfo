# Class: GRange (genomic range)

from .chromosome_range import CRange
from pyBioInfo.Utils import BlockTools


class GRange(CRange):
    def __init__(self, chrom, start=None, end=None, name=None, strand=None, blocks=None):
        if blocks:
            blocks = tuple([(block_start, block_end) for block_start, block_end in blocks])
            BlockTools.check_blocks(blocks)
            if start and start != blocks[0][0]:
                raise ValueError("Inconsistent value of start and blocks.")
            if end and end != blocks[-1][1]:
                raise ValueError("Inconsistent value of start and blocks.")
        else:
            assert start is not None and end is not None
            blocks = tuple([(start, end)])
        super(GRange, self).__init__(chrom=chrom, start=blocks[0][0], end=blocks[-1][1], name=name, strand=strand)
        self._blocks = blocks
        self._cache_length = sum([y - x for x, y in blocks])

    # Access to protected attributes.

    @property
    def blocks(self):
        return self._blocks

    # string format and length.

    def __str__(self):
        return "%s: %d-%d [strand: %s, blocks: %d]" % (
            self.__class__.__name__, 
            self._start, self._end, 
            self.strand, len(self._blocks))

    def __len__(self):
        return self._cache_length

    # index-position relevant operations.

    def get_position(self, index, strandness=True):
        length = len(self)
        if -length <= index < length:
            if strandness and self.is_reverse:
                if index >= 0:
                    index = length - 1 - index
                else:
                    index = -length - 1 - index
            return BlockTools.get_position(self.blocks, index, length=length)
        else:
            raise ValueError("Index (%d) is out of range [%d-%d]." % (index, -length, length))

    def get_index(self, position, strandness=True):
        return BlockTools.get_index(self.blocks, position, reverse=self.is_reverse and strandness)

    # Relationship

    def is_contain(self, other):
        return BlockTools.is_contain(self.blocks, other.blocks)

    def is_coincide(self, other):
        return BlockTools.is_coincide(self.blocks, other.blocks)

    def clip(self, start, end):
        blocks = BlockTools.clip(blocks=self.blocks, start=start, end=end, extend=False, check=False)
        obj = None
        if len(blocks) > 0:
            obj = GRange(chrom=self.chrom, blocks=blocks, name=self.name, strand=self.strand)
        return obj

    # formatted output.

    def format(self, fmt="bed", **kwargs):
        fmt = fmt.lower()
        line = None
        if fmt == "bed":
            ncol = kwargs.get("ncol", 12)
            cols = []
            if ncol >= 1:
                cols.append(self.chrom)
            if ncol >= 2:
                cols.append(self.start)
            if ncol >= 3:
                cols.append(self.end)
            if ncol >= 4:
                name = self.name
                if name is None:
                    name = "None"
                cols.append(name)
            if ncol >= 5:
                try:
                    score = self.score
                except AttributeError:
                    score = "."
                if score is None:
                    score = "."
                cols.append(score)
            if ncol >= 6:
                cols.append(self.strand)
            if ncol >= 7:
                cols.append(self.start)
            if ncol >= 8:
                cols.append(self.start)
            if ncol >= 9:
                try:
                    color = self.color
                except AttributeError:
                    color = "0,0,255" if self.is_forward else "255,0,0"
                cols.append(color)
            if ncol >= 10:
                cols.append(len(self.blocks))
            if ncol >= 11:
                block_sizes = [block_end - block_start for block_start, block_end in self.blocks]
                block_sizes_string = ",".join(map(str, block_sizes)) + ","
                cols.append(block_sizes_string)
            if ncol >= 12:
                block_offsets = [block_start - self.start for block_start, block_end in self.blocks]
                block_offsets_string = ",".join(map(str, block_offsets)) + ","
                cols.append(block_offsets_string)
            line = "\t".join(map(str, cols))
        elif fmt == "gff":
            cols = [self.chrom]
            try:
                source = self.source
            except AttributeError:
                source = "."
            cols.append(source)
            try:
                feature = self.feature
            except AttributeError:
                feature = "."
            cols.append(feature)

            cols.append(self.start + 1)
            cols.append(self.end)
            try:
                score = self.score
            except AttributeError:
                score = "."
            cols.append(score)

            cols.append(self.strand)
            try:
                frame = self.frame
            except AttributeError:
                frame = "."
            cols.append(frame)
            items = kwargs.get("attris")
            if items is None:
                items = ["name"]
            else:
                items = [item for item in items]
            if "name" not in items:
                items.insert(0, "name")
            array = []
            for key in items:
                if key == "name":
                    value = self.name
                else:
                    value = self[key]
                array.append("%s=%s" % (key, str(value)))
            cols.append("; ".join(array))
            line = "\t".join(map(str, cols))
        elif fmt == "tbl":
            header = kwargs.get("header")
            if header is None:
                raise RuntimeError()
            cols = []
            for key in header:
                if key == "blocks":
                    temps = []
                    for x, y in self.blocks:
                        temps.append("%d-%d" % (x, y))
                    cols.append(";".join(temps))
                else:
                    try:
                        value = self.__getattribute__(key)
                    except AttributeError:
                        value = self[key]
                    cols.append(value)
            line = "\t".join(map(str, cols))
        else:
            raise RuntimeError("Unsupported format %s" % fmt)
        return line