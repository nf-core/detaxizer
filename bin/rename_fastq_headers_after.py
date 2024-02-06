#!/usr/bin/env python
from Bio import SeqIO, bgzf
import gzip
import sys
import json
import argparse

def parse_args(args=None):

    parser = argparse.ArgumentParser(description="Renaming of the headers after running the core processes of the detaxizer pipeline.")

    parser.add_argument("-i", "--input", help="Path to input files.", type=str, required=True, nargs="+")

    parser.add_argument("-j", "--jsondict", help="Path to the gzipped json containing the IDs of the headers and the original headers.", type=str, required=True)

    parser.add_argument("-o", "--output", help="Prefix of output files.", type=str, required=True)

    return parser.parse_args()

def main():
    args = parse_args()

    fastq = args.input

    with gzip.open(args.jsondict, 'rt', encoding='utf-8') as file:
        textDict = file.read()
        headerDict = json.loads(textDict)

    if len(fastq) == 2:
        with gzip.open(fastq[0], "rt") as handle1, gzip.open(fastq[1], "rt") as handle2:
            with bgzf.BgzfWriter(args.output + "_R1_filtered.fastq.gz", "wb") as outgz1, bgzf.BgzfWriter(args.output + "_R2_filtered.fastq.gz", "wb") as outgz2:
                for i,j in zip(SeqIO.parse(handle1, "fastq"),SeqIO.parse(handle2, "fastq")):
                    id_fw = i.id
                    id_rv = j.id
                    if id_fw != id_rv:
                        sys.exit("The IDs did not match. The provided fastq files are either not sorted or a sequence is missing.")
                    lookupHeaders = headerDict[id_fw]
                    i.id = lookupHeaders[0]
                    i.description = ""
                    j.id = lookupHeaders[1]
                    j.description = ""
                    SeqIO.write(sequences=i, handle=outgz1, format="fastq")
                    SeqIO.write(sequences=j, handle=outgz2, format="fastq")
    else:
        with gzip.open(fastq[0], "rt") as handle1:
            with bgzf.BgzfWriter(args.output + "_filtered.fastq.gz", "wb") as outgz1:
                for i in SeqIO.parse(handle1, "fastq"):
                    id_fw = i.id
                    lookupHeaders = headerDict[id_fw]
                    i.id = lookupHeaders[0]
                    i.description = ""
                    SeqIO.write(sequences=i, handle=outgz1, format="fastq")

if __name__ == "__main__":
    sys.exit(main())
