import os
import gzip
import pysam
from pyBioInfo.Range import GRange
from .file import BaseFile

class VcfRecord(GRange):
    MAPPER = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "N": "N",
    }

    def __init__(self, chrom, pos, name, ref, alts, qual, filt, info):
        start = pos - 1
        end = pos
        super(VcfRecord, self).__init__(chrom=chrom,
                                        start=start, 
                                        end=end, 
                                        name=name, 
                                        strand="+")
        assert len(ref) == 1
        self.ref = ref
        if isinstance(alts, str):
            alts = [alts]
        self.alts = alts
        self.qual = qual
        self.filt = filt
        self.info = info

    def induce(self, ref, alt):
        if ref == self.ref:
            return alt in self.alts
        elif self.MAPPER[ref] == self.ref:
            return self.MAPPER[alt] in self.alts
        else:
            # assert False
            print(ref, alt, self.ref, self.alts)
            return False


class VcfFile(BaseFile):
    READ_MODE_NONE = 0
    READ_MODE_NORMAL = 1
    READ_MODE_RANDOM = 2

    def __init__(self, path, mode="r", template=None, header=None, require_index=False):
        if mode == "r":
            if require_index:
                raise NotImplementedError()
                # assert path.endswith(".vcf.gz")
                # assert os.path.exists(path + ".tbi")
        elif mode == "w":
            raise NotImplementedError()
        
        self._template = template
        self._header = header
        self._require_index = require_index
        super(VcfFile, self).__init__(path, mode)
        
        self.open()
        
        # self.path = path
        # self.mode = mode
        # self.comments = None
        # self.handle = None
        # self.read_mode = self.READ_MODE_NONE
        
        # Init comments
        # if "r" in self.mode:
        #     f = None
        #     if self.path.endswith(".vcf"):
        #         self.mode = "r"
        #         f = open(self.path, self.mode)
        #     elif self.path.endswith(".vcf.gz"):
        #         self.mode = "rt"
        #         f = gzip.open(self.path, self.mode)
        #     assert f is not None
        #     self.comments = []
        #     for line in f:
        #         line = line.strip("\n")
        #         if line.startswith("#"):
        #             self.comments.append(line)
        #         else:
        #             break

        # if "w" in self.mode:
        #     if template is None:
        #         self.comments = [
        #             "##fileformat=VCFv4.2",
        #             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
        #         ]
        #     else:
        #         self.comments = [c for c in template.comments]
                
    def open(self):
        if self._handle is None:
            if self._mode == "r":
                if self._require_index:
                    raise NotImplementedError()
                else:
                    if self._path.endswith(".gz"):
                        self._handle = gzip.open(self._path, "rt")
                    else:
                        self._handle = open(self._path)
            elif self._mode == "w":
                raise NotImplementedError()
                
    def close(self):
        if self._handle:
            self._handle.close()
        self._handle = None
        
    def __iter__(self):
        for record in self.fetch():
            yield record

    def fetch(self, chrom=None, start=None, end=None):
        if chrom:
            assert self._require_index
        else:
            for line in self._handle:
                line = line.strip("\n")
                if line.startswith("#"):
                    continue
                columns = line.split("\t")
                chrom, pos, name, ref, alts, qual, filt, info = columns
                pos = int(pos)
                alts = alts.split(",")
                qual = int(qual)
                temp = dict()
                if info != ".":
                    for item in info.split(";"):            
                        if item == "":
                            continue
                        key, value = item.split("=")
                        temp[key] = value
                info = temp
                record = VcfRecord(chrom, pos=pos, name=name, ref=ref,
                                alts=alts, qual=qual, filt=filt, info=info)
                yield record  
        
        # assert "r" in self.mode
        # if self.handle is None:
        #     assert self.read_mode == self.READ_MODE_NONE
        #     if chrom is None and start is None and end is None:
        #         if self.path.endswith(".vcf.gz"):
        #             self.handle = gzip.open(self.path, self.mode)
        #         else:
        #             self.handle = open(self.path, self.mode)
        #         self.read_mode = self.READ_MODE_NORMAL
        #     else:
        #         assert self.path.endswith(".vcf.gz")
        #         assert os.path.exists(self.path + ".tbi")
        #         self.handle = pysam.TabixFile(self.path)
        #         self.read_mode = self.READ_MODE_RANDOM

        # if chrom is None and start is None and end is None:
        #     assert self.read_mode == self.READ_MODE_NORMAL
        #     lines = self.handle
        # else:
        #     assert self.read_mode == self.READ_MODE_RANDOM
        #     lines = self.handle.fetch(chrom, start, end)

        # for line in lines:
        #     line = line.strip("\n")
        #     if line.startswith("#"):
        #         continue
        #     columns = line.split("\t")
        #     chrom, pos, name, ref, alts, qual, filt, info = columns
        #     pos = int(pos)
        #     alts = alts.split(",")
        #     qual = int(qual)
        #     temp = dict()
        #     if info != ".":
        #         for item in info.split(";"):            
        #             if item == "":
        #                 continue
        #             key, value = item.split("=")
        #             temp[key] = value
        #     info = temp
        #     record = VcfRecord(chrom, pos=pos, name=name, ref=ref,
        #                        alts=alts, qual=qual, filt=filt, info=info)
        #     yield record

    def write(self, record):
        assert "w" in self.mode
        if self.handle is None:
            if self.path.endswith(".gz"):
                self.mode = "wb"
                self.handle = gzip.open(self.path, self.mode)
            else:
                self.mode = "w+"
                self.handle = open(self.path, self.mode)
            for comment in self.comments:
                self.handle.write(comment + "\n")
        columns = [record.chrom, record.start + 1, record.name,
                   record.ref, ",".join(record.alts), record.qual, record.filt]
        info = {"SB": record.strand}
        for key, value in record.info.items():
            info[key] = value
        array = []
        for key, value in info.items():
            array.append("%s=%s" % (key, value))
        columns.append(";".join(array))
        self.handle.write("\t".join(map(str, columns)) + "\n")
