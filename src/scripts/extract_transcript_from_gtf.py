#!/usr/scripts/env python

from pyBioInfo.IO import GtfTranscriptReader

path = "/home/chenzonggui/mu01_cluster/genome/hg19/gencode.v19.annotation.sorted.gtf"
opath = "/home/chenzonggui/mu01_cluster/genome/hg19/gencode.v19.annotation.rRNA.bed"
reader = GtfTranscriptReader(path)
with open(opath, "w+") as fw:
    for transcript in reader:
        if transcript.transcript.gene_type == "rRNA":
            fw.write(transcript.format("BED") + "\n")
