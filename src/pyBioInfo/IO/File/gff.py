import os
import gzip
from collections import OrderedDict, defaultdict
import pysam
from .base import BaseFile
from pyBioInfo.Range import TRange, CRange

"""
GFF format:
1   chrom
2   source
3   feature
4   start (1-base)
5   end (1-base, not included)
6   score
7   strand
8   frame
9   attributes  ID=AAA;Parent=BBB
"""

# Record


class GffRecord(CRange):
    def __init__(self, chrom, start, end, name, strand, score, frame, source, feature, attributes):
        super(GffRecord, self).__init__(chrom=chrom,
                                        start=start,
                                        end=end,
                                        name=name,
                                        strand=strand)
        self.score = score
        self.frame = frame
        self.source = source
        self.feature = feature
        self.attributes = attributes

    def __str__(self):
        return self.format(fmt="GFF")

    def format(self, fmt="GFF"):
        fmt = fmt.upper()
        if fmt == "GFF":
            attri_str = ";".join(["%s=%s" % (k, v)
                                  for k, v in self.attributes.items()])
            values = [self.chrom, self.source, self.feature,
                      self.start + 1, self.end, self.score,
                      self.strand, self.frame, attri_str]
            return "\t".join(map(str, values))
        else:
            raise NotImplementedError()
              
                
class GffFile(BaseFile):
    def __init__(self, path, mode="r", random=None):
        if mode == "r" and random is None:
            random = path.endswith(".gff.gz") and os.path.exists(path + ".tbi")
            
        if random:
            assert mode == "r"
            assert path.endswith(".gz")
            assert os.path.exists(path)
            assert os.path.exists(path + ".tbi")
        super(GffFile, self).__init__(path, mode)
        self.comments = []
        self._random = random
        self.open()

    def open(self):
        if self._handle is None:
            if self._random:
                self._handle = pysam.TabixFile(self._path)
            else:
                if self._mode == "r":
                    if self._path.endswith(".gz"):
                        self._handle = gzip.open(self._path, "rt")
                    else:
                        self._handle = open(self._path)
                else:
                    if self._path.endswith(".gz"):
                        self._handle = gzip.open(self._path, "wt")
                    else:
                        self._handle = open(self._path, "w+")
                    
    def close(self):
        if self._handle is not None:
            self._handle.close()
            self._handle = None

    def fetch(self, chrom=None, start=None, end=None):
        if self._random:
            if chrom is None:
                for c in sorted(self._handle.contigs):
                    for line in self._handle.fetch(c):
                        yield self.parse_gff_string(line)
            else:
                for line in self._handle.fetch(chrom, start, end):
                    yield self.parse_gff_string(line)
        else:
            if chrom is None:
                for line in self._handle:
                    line = line.strip("\n")
                    if line == "":
                        continue
                    if line.startswith("#"):
                        self.comments.append(line)
                    else:
                        yield self.parse_gff_string(line)
            else:
                RuntimeError("Random access is not supported! Please set random=True.")


    @classmethod
    def parse_gff_string(cls, line):
        values = line.split("\t")
        assert len(values) == 9
        chrom = values[0]
        source = values[1]
        feature = values[2]
        start = int(values[3]) - 1
        end = int(values[4])
        score = values[5]
        strand = values[6]
        frame = values[7]
        attributes_str = values[8]
        # Make sure the strand value is valid
        assert strand == "+" or strand == "-" or strand == "."
        # Make sure the frame value is valid
        assert frame in ["0", "1", "2", "."]
        # Convert the attributes string to OrderedDict
        attributes = OrderedDict()
        for item in attributes_str.split(";"):
            item = item.strip()
            if item == "":
                continue
            x = item.find("=")
            k = item[:x]
            v = item[x + 1:].strip("\"")
            attributes[k] = v
        assert "ID" in attributes
        name = attributes["ID"]
        record = GffRecord(chrom=chrom,
                            start=start,
                            end=end,
                            name=name,
                            strand=strand,
                            score=score,
                            frame=frame,
                            source=source,
                            feature=feature,
                            attributes=attributes)
        return record
    
    
class GffTranscript(TRange):
    def __init__(self, chrom, start, end, name, strand, blocks, thick, records):
        super(GffTranscript, self).__init__(chrom=chrom,
                                            start=start,
                                            end=end,
                                            name=name,
                                            strand=strand,
                                            blocks=blocks,
                                            thick=thick)
        self.records = records  # {"feature": [record1, record2]}

    def format(self, fmt="GFF"):
        fmt = fmt.upper()
        if fmt == "GFF":
            lines = []
            for key in self.records.keys():
                for record in self.records[key]:
                    lines.append(record.format("GFF"))
            return "\n".join(lines)
        else:
            return super(GffTranscript, self).format(fmt=fmt)


class GffTranscriptBuilder(object):
    """
    The structure of feature type:
    gene -> transcript -> exon and CDS
    """

    def __init__(self, records):
        self.records = records

    @classmethod
    def construct_gff_transcript(cls, records):
        assert len(records) > 0
        transcript_id = records[0].attributes["Parent"]
        transcript_records = dict()
        for record in records:
            assert record.attributes["Parent"] == transcript_id
            feature = record.feature
            if feature not in transcript_records:
                transcript_records[feature] = list()
            transcript_records[feature].append(record)
        assert "exon" in transcript_records
        exons = transcript_records["exon"]
        chrom = exons[0].chrom
        strand = exons[0].strand
        blocks = [(exon.start, exon.end) for exon in exons]
        blocks.sort()
        start, end = blocks[0][0], blocks[-1][1]
        thick = None
        if "CDS" in transcript_records:
            cdss = transcript_records["CDS"]
            blocks_cds = [(cds.start, cds.end) for cds in cdss]
            blocks_cds.sort()
            thick = (blocks_cds[0][0], blocks_cds[-1][1])
        transcript = GffTranscript(chrom=chrom,
                                   start=start,
                                   end=end,
                                   name=transcript_id,
                                   strand=strand,
                                   blocks=blocks,
                                   thick=thick,
                                   records=transcript_records)
        return transcript

    def process_one_chromosome(self, records):
        # {ID: record, ...}
        id_to_record_dict = dict()
        # {Parent ID: [Child 1, Child 2, Child 3, ...], ...}
        parent_to_childs_dict = defaultdict(list)
        # {Transcript ID 1, Transcript ID 2, Transcript ID 3, ...}
        transcript_id_set = set()
        for record in records:
            id_to_record_dict[record.attributes["ID"]] = record
            if "Parent" in record.attributes:
                parent_to_childs_dict[record.attributes["Parent"]].append(
                    record)
            if record.feature == "transcript":
                # Make sure the transcript id is unique.
                assert record.attributes["ID"] not in transcript_id_set
                transcript_id_set.add(record.attributes["ID"])
            # if record.feature == "exon" or record.feature == "CDS":
            #     transcript_id_set.add(record.attributes["Parent"])

        transcripts = []
        for transcript_id in transcript_id_set:
            transcript = self.construct_gff_transcript(
                parent_to_childs_dict[transcript_id])
            transcript_record = id_to_record_dict[transcript_id]
            assert transcript_record.feature == "transcript" or transcript_record.feature == "mRNA"
            transcript.records[transcript_record.feature] = [transcript_record]
            if "Parent" in transcript_record.attributes:
                parent_id = transcript_record.attributes["Parent"]
                if parent_id in id_to_record_dict:
                    parent = id_to_record_dict[parent_id]
                    assert parent.feature == "gene"
                    transcript.records[parent.feature] = [parent]
            transcripts.append(transcript)
        transcripts.sort()
        return transcripts

    def __iter__(self):
        chrom = None    # Current chromsome
        array = None    # Records of current chromosome
        for record in self.records:
            if chrom is None:
                chrom = record.chrom
                array = [record]
            else:
                if record.chrom == chrom:
                    array.append(record)
                elif record.chrom > chrom:
                    for transcript in self.process_one_chromosome(array):
                        yield transcript
                    chrom = record.chrom
                    array = [record]
                else:
                    raise RuntimeError("Records is not sorted by chrom!")
        if chrom:
            for transcript in self.process_one_chromosome(array):
                yield transcript
            chrom = None
            array = None

# Gene


class GffGene(TRange):
    def __init__(self, chrom, start, end, name, strand, blocks, thick, transcripts, records):
        super(GffGene, self).__init__(chrom=chrom,
                                      start=start,
                                      end=end,
                                      strand=strand,
                                      name=name,
                                      blocks=blocks,
                                      thick=thick)
        # [Transcript 1, Transcript 2, Transcript 3, ...]
        self.transcripts = transcripts
        # {
        #   "gene": [record],
        #   "transcript": [record1, record2, ...],
        #   "exon": [record1, record2, ...],
        #   "CDS": [record1, record2, ...],
        #   ....
        # }
        self.records = records


class GffGeneBuilder(object):
    def __init__(self, records):
        self.records = records

    def process_one_chromosome(self, records):
        id_to_record_dict = dict()
        parent_child_dict = defaultdict(list)
        gene_id_set = set()
        for record in records:
            id_to_record_dict[record.attributes["ID"]] = record
            if record.feature == "gene":
                gene_id_set.add(record.attributes["ID"])
            if "Parent" in record.attributes:
                parent_child_dict[record.attributes["Parent"]].append(record)

        genes = []
        for gene_id in gene_id_set:
            gene_record = id_to_record_dict[gene_id]
            assert gene_record.feature == "gene"
            gene_records = defaultdict(list)
            gene_records[gene_record.feature] = [gene_record]

            transcripts = []
            for transcript_record in parent_child_dict[gene_id]:
                assert transcript_record.feature == "transcript" or transcript_record.feature == "mRNA"
                gene_records[transcript_record.feature].append(transcript_record)
                transcript_id = transcript_record.attributes["ID"]
                transcript_record_list = parent_child_dict[transcript_id]
                for record in transcript_record_list:
                    gene_records[record.feature].append(record)
                transcript = GffTranscriptBuilder.construct_gff_transcript(transcript_record_list)
                transcript.records[transcript_record.feature] = [
                    transcript_record]
                transcript.records[gene_record.feature] = [gene_record]
                transcripts.append(transcript)
            transcripts.sort()

            chrom = transcripts[0].chrom
            start = min([transcript.start for transcript in transcripts])
            assert start == gene_record.start
            end = max([transcript.end for transcript in transcripts])
            assert end == gene_record.end
            strand = transcripts[0].strand
            blocks = None
            thick = None
            gene = GffGene(chrom=chrom,
                           start=start,
                           end=end,
                           strand=strand,
                           name=gene_id,
                           records=gene_records,
                           transcripts=transcripts,
                           blocks=blocks, thick=thick)
            genes.append(gene)
        genes.sort()
        return genes

    def __iter__(self):
        chrom = None
        array = None
        for record in self.records:
            if chrom is None:
                chrom = record.chrom
                array = [record]
            else:
                if record.chrom == chrom:
                    array.append(record)
                elif record.chrom > chrom:
                    for gene in self.process_one_chromosome(array):
                        yield gene
                    chrom = record.chrom
                    array = [record]
                else:
                    raise RuntimeError("Records is not sorted by chrom!")
        if chrom:
            for gene in self.process_one_chromosome(array):
                yield gene
            chrom = None
            array = None


class GffTreeBuilder(object):
    def __init__(self, records):
        self.records = records

    def process_one_chromosome(self, records):
        id_to_record_dict = dict()
        parent_child_dict = defaultdict(list)
        no_parent_records = []
        for record in records:
            assert "ID" in record.attributes.keys()
            id_to_record_dict[record.attributes["ID"]] = record
            if "Parent" in record.attributes:
                parent_child_dict[record.attributes["Parent"]].append(record)
            else:
                no_parent_records.append(record)
        for pid, items in parent_child_dict.items():
            parent = id_to_record_dict[pid]
            parent.childs = parent_child_dict[pid]
        
        array = []
        for record in no_parent_records:
            array.append(record)
        array.sort()
        return array

    def __iter__(self):
        chrom = None
        array = None
        for record in self.records:
            if chrom is None:
                chrom = record.chrom
                array = [record]
            else:
                if record.chrom == chrom:
                    array.append(record)
                elif record.chrom > chrom:
                    for gene in self.process_one_chromosome(array):
                        yield gene
                    chrom = record.chrom
                    array = [record]
                else:
                    print(record.chrom, chrom, record.chrom < chrom)
                    raise RuntimeError("Records is not sorted by chrom!")
        if chrom:
            for gene in self.process_one_chromosome(array):
                yield gene
            chrom = None
            array = None