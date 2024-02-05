#!/usr/bin/env python
from Bio import SeqIO, bgzf
import gzip
import sys
import json
import argparse

def parse_args(args=None):

    parser = argparse.ArgumentParser(description="Renaming of the headers before running the core processes of the detaxizer pipeline.")

    parser.add_argument("-i", "--input", help="Path to input files.", type=str, required=True, nargs="+")

    parser.add_argument("-o", "--output", help="Prefix of output files.", type=str, required=True)

    return parser.parse_args()

def renameReadsPaired(reads):
    read_fw = reads[0]
    read_rv = reads[1]
    read_dict = {}
    if "/1" in read_fw and "/2" in read_rv and " " not in read_fw and " " not in read_rv:
        read_fw_stripped = read_fw.strip("1").strip("/")
        read_rv_stripped = read_rv.strip("2").strip("/")
        if read_fw_stripped != read_rv_stripped:
            sys.exit("Read IDs were not matching! Please provide matching headers.")
        else:
            read_dict[read_fw_stripped] = [read_fw, read_rv]
            read_fw_stripped = read_fw_stripped + " 1:N:10:"
            read_rv_stripped = read_rv_stripped + " 2:N:10:"
            read_renamed = [read_fw_stripped,read_rv_stripped]
    elif "/1" in read_fw and "/2" in read_rv and " " in read_fw and " " in read_rv:
        read_fw_stripped = read_fw.strip("1").strip("/").split(" ")[0]
        read_rv_stripped = read_rv.strip("2").strip("/").split(" ")[0]
        if read_fw_stripped != read_rv_stripped:
            sys.exit("Read IDs were not matching! Please provide matching headers.")
        else:
            read_dict[read_fw_stripped] = [read_fw, read_rv]
            read_fw_stripped = read_fw_stripped + " 1:N:10:"
            read_rv_stripped = read_rv_stripped + " 2:N:10:"
            read_renamed = [read_fw_stripped,read_rv_stripped]
    elif " 1:" in read_fw and " 2:" in read_rv:
        read_fw_stripped = read_fw.split(" ")[0]
        read_rv_stripped = read_rv.split(" ")[0]
        if read_fw_stripped != read_rv_stripped:
            sys.exit("Read IDs were not matching! Please provide matching headers.")
        else:
            read_dict[read_fw_stripped] = [read_fw, read_rv]
            read_fw_stripped = read_fw_stripped + " 1:N:10:"
            read_rv_stripped = read_rv_stripped + " 2:N:10:"
            read_renamed = [read_fw_stripped,read_rv_stripped]
    else:
        sys.exit("The headers were not matching the patterns!")
    return (read_dict,read_renamed)

def renameReadSingle(read):
    read_dict = {}
    if "/1" in read and " " not in read:
        read_fw_stripped = read.strip("1").strip("/")
        read_dict[read_fw_stripped] = [read]
        read_fw_stripped = read_fw_stripped + " 1:N:10:"
        read_renamed = [read_fw_stripped]
    elif "/1" in read and " " in read:
        read_fw_stripped = read.strip("1").strip("/").split(" ")[0]
        read_dict[read_fw_stripped] = [read]
        read_fw_stripped = read_fw_stripped + " 1:N:10:"
        read_renamed = [read_fw_stripped]
    elif " 1:" in read:
        read_fw_stripped = read.split(" ")[0]
        read_dict[read_fw_stripped] = [read]
        read_fw_stripped = read_fw_stripped + " 1:N:10:"
        read_renamed = [read_fw_stripped]
    elif " " in read:
        read_fw_stripped = read.split(" ")[0]
        read_dict[read_fw_stripped] = [read]
        read_fw_stripped = read_fw_stripped + " 1:N:10:"
        read_renamed = [read_fw_stripped]
    else:
        sys.exit("The headers were not matching the patterns!")
    return (read_dict,read_renamed)

def main():
    args = parse_args()

    fastq = args.input

    if len(fastq) == 2:
        renamed = {}
        with gzip.open(fastq[0], "rt") as handle1, gzip.open(fastq[1], "rt") as handle2:
            with bgzf.BgzfWriter(args.output + "_R1_renamed.fastq.gz", "wb") as outgz1, bgzf.BgzfWriter(args.output + "_R2_renamed.fastq.gz", "wb") as outgz2:
                for i,j in zip(SeqIO.parse(handle1, "fastq"),SeqIO.parse(handle2, "fastq")):
                    header_fw = i.description
                    header_rv = j.description
                    headers = (header_fw,header_rv)
                    headers = renameReadsPaired(headers)
                    renamed.update(headers[0])
                    i.description = headers[1][0]
                    i.id = headers[1][0].split(" ")[0]
                    j.description = headers[1][1]
                    j.id = headers[1][1].split(" ")[0]
                    SeqIO.write(sequences=i, handle=outgz1, format="fastq")
                    SeqIO.write(sequences=j, handle=outgz2, format="fastq")
        with gzip.open(args.output + "_headers.json.gz", "wt", encoding="utf-8") as outfile:
            json.dump(renamed, outfile)
    else:
        renamed = {}
        with gzip.open(fastq[0], "rt") as handle1:
            with bgzf.BgzfWriter(args.output + "_renamed.fastq.gz", "wb") as outgz1:
                for i in SeqIO.parse(handle1, "fastq"):
                    header_fw = i.description
                    headers = renameReadSingle(header_fw)
                    renamed.update(headers[0])
                    i.description = headers[1][0]
                    i.id = headers[1][0].split(" ")[0]
                    SeqIO.write(sequences=i, handle=outgz1, format="fastq")
        with gzip.open(args.output + "_headers.json.gz", "wt", encoding="utf-8") as outfile:
            json.dump(renamed, outfile)

if __name__ == "__main__":
    sys.exit(main())
