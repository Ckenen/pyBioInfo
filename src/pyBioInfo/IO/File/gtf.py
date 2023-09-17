from collections import defaultdict, OrderedDict
import gzip
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


class GtfFile(object):
    def __init__(self, path, mode="r"):
        self.path = path
        self.mode = mode
        self.comments = []
        self.handle = None

    def open(self):
        assert self.handle is None
        if self.path.endswith(".gz"):
            if self.mode == "r":
                self.mode = "rt"
            self.handle = gzip.open(self.path, self.mode)
        else:
            self.handle = open(self.path, self.mode)

    def close(self):
        if self.handle:
            self.handle.close()
            self.handle = None

    def write(self):
        raise NotImplementedError()

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __iter__(self):
        for line in self.handle:
            line = line.strip("\n")
            if line.startswith("#"):
                self.comments.append(line)
                continue
            values = line.split("\t")
            # print(line)
            # print(values)
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
            # if frame in ["0", "1", "2"]:
            #     frame = int(frame)
            # elif frame == ".":
            #     frame = None
            # else:
            #     raise ValueError()

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
            yield record

# Transcript


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
                continue
            container[record.attributes["transcript_id"]].append(record)

        transcripts = []
        for transcript_id in container.keys():
            transcript = self.construct_gtf_transcript(
                container[transcript_id])
            gene_id = transcript.records["transcript"][0].attributes["gene_id"]
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

# Gene


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
