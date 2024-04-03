import os
from collections import defaultdict, OrderedDict
import gzip
import pysam
from .file import BaseFile
from pyBioInfo.Range import GRange, CRange, TRange
from pyBioInfo.Utils import BlockTools

"""
GTF format:
colum-number    content    values/format
1	chromosome name	chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,M} or GRC accession a
2	annotation source	{ENSEMBL,HAVANA}
3	feature type	{gene,transcript,exon,CDS,UTR,start_codon,stop_codon,Selenocysteine}
4	genomic start location	integer-value (1-based)
5	genomic end location	integer-value
6	score(not used)	.
7	genomic strand	{+,-}
8	genomic phase (for CDS features)	{0,1,2,.}
9	additional information as key-value pairs	see below (GENCODE)

reference: https://www.gencodegenes.org/pages/data_format.html
"""

# Record


class GtfRecord(CRange):
    def __init__(self, chrom, start, end, strand, score, frame, source, feature, attributes):
        super(GtfRecord, self).__init__(chrom=chrom,
                                        start=start,
                                        end=end,
                                        strand=strand)
        self.score = score
        self.frame = frame
        self.source = source
        self.feature = feature
        self.attributes = attributes

    def format(self, fmt="GTF"):
        fmt = fmt.upper()
        if fmt == "GTF":
            attri_str = "; ".join(["%s \"%s\"" % (k, v)
                                   for k, v in self.attributes.items()]) + ";"
            values = [self.chrom, self.source, self.feature,
                      self.start + 1, self.end, self.score,
                      self.strand, self.frame, attri_str]
            return "\t".join(map(str, values))
        else:
            raise NotImplementedError()
        

class GtfFile(BaseFile):
    def __init__(self, path, mode="rt"):
        super(GtfFile, self).__init__(path, mode)
        self.comments = []
        self.open()

    def open(self):
        if self._handle is None:
            if self._mode == "rt":
                if self._path.endswith(".gz"):
                    self._handle = gzip.open(self._path, "rt")
                else:
                    self._handle = open(self._path, "r")
            else:
                RuntimeError("Random access is not supported! Please try GtfFileRandom.")

    def close(self):
        if self._handle is not None:
            self._handle.close()
            self._handle = None

    def fetch(self, chrom=None, start=None, end=None):
        if chrom is None:
            for line in self._handle:
                line = line.strip("\n")
                if line == "":
                    continue
                if line.startswith("#"):
                    self.comments.append(line)
                else:
                    yield self.parse_gtf_string(line)
        else:
            raise RuntimeError()

    def __iter__(self):
        for x in self.fetch():
            yield x

    @classmethod
    def parse_gtf_string(cls, line):
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
        assert strand == "+" or strand == "-" or strand == "."
        assert frame in ["0", "1", "2", "."]
        attributes = OrderedDict()
        for item in attributes_str.split(";"):
            item = item.strip()
            if item == "":
                continue
            x = item.find(" ")
            k = item[:x]
            v = item[x + 1:].strip("\"")
            attributes[k] = v
        assert "gene_id" in attributes
        if feature != "gene":
            assert "transcript_id" in attributes
        record = GtfRecord(chrom=chrom,
                            start=start,
                            end=end,
                            strand=strand,
                            score=score,
                            frame=frame,
                            source=source,
                            feature=feature,
                            attributes=attributes)
        return record


class GtfFileRandom(GtfFile):
    def __init__(self, path):
        assert path.endswith(".gz")
        assert os.path.exists(path)
        assert os.path.exists(path + ".tbi")
        super(GtfFileRandom, self).__init__(path, "rt")

    def open(self):
        assert self._mode == "rt"
        if self._handle is None:
            self._handle = pysam.TabixFile(self._path)

    def fetch(self, chrom=None, start=None, end=None):
        if chrom is None:
            for c in sorted(self._handle.contigs):
                for line in self._handle.fetch(c):
                    yield GtfFile.parse_gtf_string(line)
        else:
            for line in self._handle.fetch(chrom, start, end):
                yield GtfFile.parse_gtf_string(line)


class GtfTranscript(TRange):
    def __init__(self, chrom, start, end, name, strand, blocks, thick, records):
        super(GtfTranscript, self).__init__(chrom=chrom,
                                            start=start,
                                            end=end,
                                            name=name,
                                            strand=strand,
                                            blocks=blocks,
                                            thick=thick)
        self.records = records  # {"feature": [record1, record2]}


class GtfTranscriptBuilder(object):
    def __init__(self, records):
        self.records = records

    @classmethod
    def construct_gtf_transcript(cls, records):
        assert len(records) > 0
        transcript_id = records[0].attributes["transcript_id"]
        transcript_records = defaultdict(list)
        for record in records:
            assert record.attributes["transcript_id"] == transcript_id
            transcript_records[record.feature].append(record)
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
        # print(transcript_id)
        # print(blocks)
        transcript = GtfTranscript(chrom=chrom,
                                   start=start,
                                   end=end,
                                   name=transcript_id,
                                   strand=strand,
                                   blocks=blocks,
                                   thick=thick,
                                   records=transcript_records)
        return transcript

    def process_one_chromosome(self, records):
        # {"gene_id": record, "gene_id": record, ...}
        genes = dict()
        # {"transcript_id": [record 1, record 2, record 3, ...]}
        container = defaultdict(list)
        for record in records:
            if record.feature == "gene":
                genes[record.attributes["gene_id"]] = record
            else:
                container[record.attributes["transcript_id"]].append(record)

        transcripts = []
        for transcript_id in container.keys():
            transcript = self.construct_gtf_transcript(
                container[transcript_id])
            gene_id = container[transcript_id][0].attributes["gene_id"]
            if gene_id in genes.keys():
                assert "gene" not in transcript.records.keys()
                transcript.records["gene"] = [genes[gene_id]]
            transcripts.append(transcript)
        transcripts.sort()
        return transcripts

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
                else:
                    for transcript in self.process_one_chromosome(array):
                        yield transcript
                    chrom = record.chrom
                    array = [record]
        if chrom:
            for transcript in self.process_one_chromosome(array):
                yield transcript
            chrom = None
            array = None


class GtfGene(TRange):
    def __init__(self, chrom, start, end, name, strand, blocks, thick, transcripts, records):
        super(GtfGene, self).__init__(chrom=chrom, start=start, end=end,
                                      strand=strand, name=name, blocks=blocks, thick=thick)
        # [Transcript 1, Transcript 2, Transcript 3, ...]
        self.transcripts = transcripts
        # {
        #   "gene": [record],
        #   "transcript": [record 1, record 2, record 3, ...],
        #   "exon": [record 1, record 2, record 3, ...]
        # }
        self.records = records


class GtfGeneBuilder(object):
    def __init__(self, records):
        self.records = records

    def process_one_chromosome(self, records):
        # {"gene_id": [record1, record2, record3]}
        container = defaultdict(list)
        for record in records:
            container[record.attributes["gene_id"]].append(record)

        genes = []
        for gene_id in container.keys():
            array = container[gene_id]
            # {
            #   "gene": [record],
            #   "transcript": [record 1, record 2, ...],
            #   "exon": [record 1, record 2, ...],
            #   "CDS": [record 1, record 2, ...],
            #   ...
            # }
            gene_records = defaultdict(list)
            for record in array:
                gene_records[record.feature].append(record)

            # {"transcript_id": [record1, record2, record3, ...]}
            container1 = defaultdict(list)
            transcripts = []
            for record in array:
                if record.feature == "gene":
                    continue
                container1[record.attributes["transcript_id"]].append(record)

            for transcript_id in container1:
                transcript = GtfTranscriptBuilder.construct_gtf_transcript(
                    container1[transcript_id])
                if "gene" in gene_records:
                    transcript.records["gene"] = [gene_records["gene"]]
                transcripts.append(transcript)
            transcripts.sort()
            chrom = transcripts[0].chrom
            start = min([transcript.start for transcript in transcripts])
            end = max([transcript.end for transcript in transcripts])
            strand = transcripts[0].strand
            blocks = []
            for transcript in transcripts:
                for block in transcript.blocks:
                    blocks.append(block)
            blocks = list(sorted(blocks))
            blocks = BlockTools.suture(blocks)
            
            thick = None
            gene = GtfGene(chrom=chrom,
                           start=start,
                           end=end,
                           strand=strand,
                           name=gene_id,
                           records=gene_records,
                           transcripts=transcripts,
                           blocks=blocks,
                           thick=thick)
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
                else:
                    for gene in self.process_one_chromosome(array):
                        yield gene
                    chrom = record.chrom
                    array = [record]
        if chrom:
            for gene in self.process_one_chromosome(array):
                yield gene
            chrom = None
            array = None
