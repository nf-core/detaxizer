#!/usr/bin/env python

# Written by Jannik Seidel and released under MIT license.

from Bio import SeqIO, bgzf
import gzip
import sys
import argparse
import re

def parse_args(args=None):

    parser = argparse.ArgumentParser(description="Renaming of the headers before running the core processes of the detaxizer pipeline.")

    parser.add_argument("-i", "--input", help="Path to input files.", type=str, required=True, nargs="+")

    parser.add_argument("-o", "--output", help="Prefix of output files.", type=str, required=True)

    return parser.parse_args()

def renameReadsPaired(reads: tuple, filenames: str) -> tuple:
    read_fw = reads[0]
    read_rv = reads[1]
    read_dict = {}
    # Illumina format pre-CASAVA 1.8
    pattern1 = r"^\S+/[1,2]$"       # matches for example the following header 'example.1/1'
    # Illumina format pre-CASAVA 1.8 with additional information
    pattern2 = r"^\S+/[1,2]\s"      # matches for example the following header 'example.1/1 additionalInformation'
    # Illumina format post-CASAVA 1.8
    pattern3 = r"^\S+\s\S+$"        # matches for example the following header 'readID1 additionalTechnicalInformation'
    # Illumina format post-CASAVA 1.8 with additional Information
    pattern4 = r"^\S+\s\S+\s\S+"    # matches for example the following header 'readID1 additionalTechnicalInformation additionalInformation'
    # Any other pattern without spaces
    pattern5 = r"^\S+$"             # matches for example the following header 'readID1'
    if bool(re.match(pattern1,read_fw)) and bool(re.match(pattern1,read_rv)):
        if read_fw.endswith("/1"):
            read_fw_stripped = read_fw[:-2]
        else:
            raise ValueError("Please provide the forward reads in short_reads_fastq_1 (where the headers are as follows: 'example.1/1').")

        if read_rv.endswith("/2"):
            read_rv_stripped = read_rv[:-2]
        else:
            raise ValueError("Please provide the reverse reads in short_reads_fastq_2 (where the headers are as follows: 'example.1/2').")

        if read_fw_stripped != read_rv_stripped:
            msg = f"Read IDs were not matching! Please provide matching IDs in the headers. The problematic reads were {read_fw} and {read_rv} in the files {filenames}."
            raise ValueError(msg)
        else:
            read_dict[read_fw_stripped] = [read_fw, read_rv]
            read_renamed = [read_fw_stripped,read_rv_stripped]
    elif bool(re.match(pattern2,read_fw)) and bool(re.match(pattern2,read_rv)):
        read_fw_split = read_fw.split(" ")[0]
        read_rv_split = read_rv.split(" ")[0]
        if read_fw_split.endswith("/1"):
            read_fw_stripped = read_fw_split[:-2]
        else:
            raise ValueError("Please provide the forward reads in short_reads_fastq_1 (where the headers are as follows: 'example.1/1 additionalInformation').")

        if read_rv_split.endswith("/2"):
            read_rv_stripped = read_rv_split[:-2]
        else:
            raise ValueError("Please provide the reverse reads in short_reads_fastq_2 (where the headers are as follows: 'example.1/2 additionalInformation').")

        if read_fw_stripped != read_rv_stripped:
            msg = f"Read IDs were not matching! Please provide matching IDs in the headers. The problematic reads were {read_fw} and {read_rv} in the files {filenames}."
            raise ValueError(msg)
        else:
            read_dict[read_fw_stripped] = [read_fw, read_rv]
            read_renamed = [read_fw_stripped,read_rv_stripped]
    elif bool(re.match(pattern3,read_fw)) and bool(re.match(pattern3,read_rv)):
        read_fw_split = read_fw.split(" ")[0]
        read_rv_split = read_rv.split(" ")[0]

        if read_fw_split != read_rv_split:
            msg = f"Read IDs were not matching! Please provide matching IDs in the headers. The problematic reads were {read_fw} and {read_rv} in the files {filenames}."
            raise ValueError(msg)
        else:
            read_dict[read_fw_split] = [read_fw, read_rv]
            read_renamed = [read_fw_split,read_rv_split]
    elif bool(re.match(pattern4,read_fw)) and bool(re.match(pattern4,read_rv)):
        read_fw_split = read_fw.split(" ")[0]
        read_rv_split = read_rv.split(" ")[0]

        if read_fw_split != read_rv_split:
            msg = f"Read IDs were not matching! Please provide matching IDs in the headers. The problematic reads were {read_fw} and {read_rv} in the files {filenames}."
            raise ValueError(msg)
        else:
            read_dict[read_fw_split] = [read_fw, read_rv]
            read_renamed = [read_fw_split,read_rv_split]
    elif bool(re.match(pattern5,read_fw)) and bool(re.match(pattern5,read_rv)):
        if read_fw != read_rv:
            msg = f"Read IDs were not matching! Please provide matching IDs in the headers. The problematic reads were {read_fw} and {read_rv} in the files {filenames}."
            raise ValueError(msg)
        else:
            read_dict[read_fw] = [read_fw, read_rv]
            read_renamed = [read_fw,read_rv]
    else:
        msg = f"The provided files, {filenames}, contained reads with headers not supported by the pipeline.\n  Please use one of the formats:\n    example.1/1\n    example.1/1 additionalInformation\n    readID1 additionalTechnicalInformation\n    readID1 additionalTechnicalInformation additionalInformation\n    readID1\nAny other format is not supported."
        raise ValueError(msg)
    return (read_dict,read_renamed)

def renameReadSingle(read: str, filename: str) -> tuple:
    read_dict = {}
    # Illumina format pre-CASAVA 1.8
    pattern1 = r"^\S+/[1,2]$"       # matches for example the following header 'example.1/1'
    # Illumina format pre-CASAVA 1.8 with additional information
    pattern2 = r"^\S+/[1,2]\s"      # matches for example the following header 'example.1/1 additionalInformation'
    # Illumina format post-CASAVA 1.8
    pattern3 = r"^\S+\s\S+$"        # matches for example the following header 'readID1 additionalTechnicalInformation'
    # Illumina format post-CASAVA 1.8 with additional Information
    pattern4 = r"^\S+\s\S+\s\S+"    # matches for example the following header 'readID1 additionalTechnicalInformation additionalInformation'
    # Any other pattern without spaces
    pattern5 = r"^\S+$"             # matches for example the following header 'readID1'
    if bool(re.match(pattern1,read)):
        read_stripped = read[:-2]
        read_dict[read_stripped] = [read]
        read_renamed = [read_stripped]
    elif bool(re.match(pattern2,read)):
        read_split = read.split(" ")[0]
        read_stripped = read_split[:-2]
        read_dict[read_stripped] = [read]
        read_renamed = [read_stripped]
    elif bool(re.match(pattern3,read)):
        read_split = read.split(" ")[0]
        read_dict[read_split] = [read]
        read_renamed = [read_split]
    elif bool(re.match(pattern4,read)):
        read_split = read.split(" ")[0]
        read_dict[read_split] = [read]
        read_renamed = [read_split]
    elif bool(re.match(pattern5,read)):
            read_dict[read] = [read]
            read_renamed = [read]
    else:
        msg = f"The provided file, {filename}, contained reads with headers not supported by the pipeline.\n  Please use one of the formats:\n    example.1/1\n    example.1/1 additionalInformation\n    readID1 additionalTechnicalInformation\n    readID1 additionalTechnicalInformation additionalInformation\n    readID1\nAny other format is not supported."
        raise ValueError(msg)
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
                    renamedHeaders = renameReadsPaired(headers, " and ".join(fastq))
                    renamed.update(renamedHeaders[0])
                    i.description = renamedHeaders[1][0]
                    i.id = renamedHeaders[1][0]
                    j.description = renamedHeaders[1][1]
                    j.id = renamedHeaders[1][1]
                    SeqIO.write(sequences=i, handle=outgz1, format="fastq")
                    SeqIO.write(sequences=j, handle=outgz2, format="fastq")
        with gzip.open(args.output + "_headers_fw.txt.gz", "wt", encoding="utf-8") as outfile1, gzip.open(args.output + "_headers_rv.txt.gz", "wt", encoding="utf-8") as outfile2:
            for key in renamed.keys():
                outfile1.write(key + "\t" + renamed[key][0] + "\n")
                outfile2.write(key + "\t" + renamed[key][1] + "\n")
    else:
        renamed = {}
        with gzip.open(fastq[0], "rt") as handle1:
            with bgzf.BgzfWriter(args.output + "_renamed.fastq.gz", "wb") as outgz1:
                for i in SeqIO.parse(handle1, "fastq"):
                    header_fw = i.description
                    renamedHeader = renameReadSingle(header_fw, fastq)
                    renamed.update(renamedHeader[0])
                    i.description = renamedHeader[1][0]
                    i.id = renamedHeader[1][0]
                    SeqIO.write(sequences=i, handle=outgz1, format="fastq")
        with gzip.open(args.output + "_headers.txt.gz", "wt", encoding="utf-8") as outfile:
            for key in renamed.keys():
                outfile.write(key + "\t" + renamed[key][0] + "\n")

if __name__ == "__main__":
    sys.exit(main())
