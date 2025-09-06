#!/usr/bin/env python
import sys, os, subprocess, optparse
from collections import defaultdict
import pandas as pd
from pyBioInfo.IO.File import GtfFile, GtfTranscriptBuilder


def output_gene_info(records, outfile):
    rows = [] # gene
    for record in records:
        if record.feature == "gene":
            row = [
                record.attributes["gene_id"],
                record.attributes["gene_type"],
                record.attributes["gene_name"],
                record.chrom,
                record.start,
                record.end,
                record.strand,
            ]
            rows.append(row)
    columns = ["GeneID", "GeneType", "GeneName", "Chrom", "Start", "End", "Strand"]
    m = pd.DataFrame(rows, columns=columns)
    m.to_csv(outfile, sep="\t" if outfile.endswith(".tsv") else ",", index=False)


def output_transcript_info(records, outfile):
    lengths = defaultdict(int)
    rows = [] # transcript
    for record in records:
        if record.feature == "transcript":
            row = [
                record.attributes["transcript_id"],
                record.attributes["transcript_type"],
                record.attributes["transcript_name"],
                record.attributes["gene_id"],
                record.attributes["gene_type"],
                record.attributes["gene_name"],
                record.chrom,
                record.start,
                record.end,
                record.strand,
            ]
            rows.append(row)
        elif record.feature == "exon":
            lengths[record.attributes["transcript_id"]] += len(record)
    columns = ["TranscriptID", "TranscriptType", "TranscriptName", "GeneID", "GeneType", "GeneName", "Chrom", "Start", "End", "Strand"]
    m = pd.DataFrame(rows, columns=columns)
    m["Length"] = m["TranscriptID"].map(lengths)
    m.to_csv(outfile, sep="\t" if outfile.endswith(".tsv") else ",", index=False)


def output_gene_bed(records, outfile, fullname=False):
    assert outfile.endswith(".bed.gz")
    tmpfile = outfile[:-3]
    with open(tmpfile, "w+") as fw:
        for record in records:
            if record.feature == "gene":
                row = [
                    record.chrom,
                    record.start,
                    record.end,
                    ":".join([
                        record.attributes["gene_id"], 
                        record.attributes["gene_name"],
                        record.attributes["gene_type"],
                    ]) if fullname else record.attributes["gene_id"],
                    ".",
                    record.strand,
                ]
                fw.write("\t".join(map(str, row)) + "\n")
    cmd1 = "bgzip -o %s %s" % (outfile, tmpfile)
    cmd2 = "tabix -p bed %s" % outfile
    subprocess.check_call(cmd1, shell=True)
    subprocess.check_call(cmd2, shell=True)
    os.remove(tmpfile)


def output_transcript_bed(records, outfile, fullname=False):
    assert outfile.endswith(".bed.gz")
    tmpfile = outfile[:-3]
    transcripts = GtfTranscriptBuilder(records)
    transcripts = list(sorted(transcripts, key=lambda t: (t.chrom, t.start, t.end, t.name)))
    with open(tmpfile, "w+") as fw:
        for transcript in transcripts:
            record = transcript.records["transcript"][0]
            if fullname:
                transcript.name = ":".join([
                    record.attributes["gene_id"],
                    record.attributes["gene_name"],
                    record.attributes["gene_type"],
                    record.attributes["transcript_id"],
                    record.attributes["transcript_name"],
                    record.attributes["transcript_type"],
                ])
            fw.write(transcript.format("bed") + "\n")
    cmd1 = "bgzip -o %s %s" % (outfile, tmpfile)
    cmd2 = "tabix -p bed %s" % outfile
    subprocess.check_call(cmd1, shell=True)
    subprocess.check_call(cmd2, shell=True)
    os.remove(tmpfile)
    
    
def load_records(gtf):
    records = []
    features = ["gene", "transcript", "exon", "CDS"]
    with GtfFile(gtf) as f:
        for record in f:
            if record.feature not in features:
                continue
            records.append(record)
    return records
        
        
def main():
    parser = optparse.OptionParser(usage="""
                                   
    %prog input.gtf prefix

This script will output following files:
    1. prefix.gene_info.csv
    2. prefix.transcript_info.csv
    3. prefix.genes.bed.gz
    4. prefix.genes.fullname.bed.gz
    5. prefix.transcripts.bed.gz
    6. prefix.transcripts.fullname.bed.gz""")
    
    _, args = parser.parse_args()
    
    gtf, prefix = args
        
    gene_info = prefix + ".gene_info.csv"
    transcript_info = prefix + ".transcript_info.csv"
    gene_bed = prefix + ".genes.bed.gz"
    transcript_bed = prefix + ".transcripts.bed.gz"
    gene_full_name_bed = prefix + ".genes.fullname.bed.gz"
    transcript_full_name_bed = prefix + ".transcripts.fullname.bed.gz"
    
    records = load_records(gtf)
    
    output_gene_info(records, gene_info)
    output_transcript_info(records, transcript_info)
    output_gene_bed(records, gene_bed)
    output_transcript_bed(records, transcript_bed)
    output_gene_bed(records, gene_full_name_bed, fullname=True)
    output_transcript_bed(records, transcript_full_name_bed, fullname=True)    


if __name__ == "__main__":
    main()