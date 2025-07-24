#!/usr/bin/env python

# Written by Jannik Seidel and released under MIT license.
# Modified for low memory buffered writing

import dnaio
import gzip
import sys
import argparse
import re
from io import StringIO

# Number of reads to buffer before writing
BUFFER_SIZE = 100000

def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Renaming of the headers before running the core processes of the detaxizer pipeline.")
    parser.add_argument("-i", "--input", help="Path to input files.", type=str, required=True, nargs="+")
    parser.add_argument("-o", "--output", help="Prefix of output files.", type=str, required=True)
    parser.add_argument("-b", "--buffer", help="Number of reads to buffer before writing to disk.", type=int, default=BUFFER_SIZE)
    return parser.parse_args()

def process_read_header(header, is_forward=None, filename=None):
    """Process a single read header and return the renamed version."""
    # Illumina format pre-CASAVA 1.8
    pattern1 = r"^\S+/[1,2]$"
    # Illumina format pre-CASAVA 1.8 with additional information
    pattern2 = r"^\S+/[1,2]\s"
    # Illumina format post-CASAVA 1.8
    pattern3 = r"^\S+\s\S+$"
    # Illumina format post-CASAVA 1.8 with additional Information
    pattern4 = r"^\S+\s\S+\s\S+"
    # Any other pattern without spaces
    pattern5 = r"^\S+$"
    
    if bool(re.match(pattern1, header)):
        return header[:-2]  # Remove /1 or /2
    elif bool(re.match(pattern2, header)):
        return header.split(" ")[0][:-2]  # Remove /1 or /2 from first part
    elif bool(re.match(pattern3, header)) or bool(re.match(pattern4, header)):
        return header.split(" ")[0]  # Return first part before space
    elif bool(re.match(pattern5, header)):
        return header  # Return as is
    else:
        msg = f"The provided file contained reads with headers not supported by the pipeline.\n  Please use one of the formats:\n    example.1/1\n    example.1/1 additionalInformation\n    readID1 additionalTechnicalInformation\n    readID1 additionalTechnicalInformation additionalInformation\n    readID1\nAny other format is not supported."
        raise ValueError(msg)

def write_buffer(buffer, output_file):
    """Write buffer contents to file and clear the buffer."""
    if buffer.getvalue():
        output_file.write(buffer.getvalue())
        buffer.seek(0)
        buffer.truncate(0)

def main():
    args = parse_args()
    fastq = args.input
    buffer_size = args.buffer
    
    if len(fastq) == 2:
        # Create buffers for header mappings
        fw_buffer = StringIO()
        rv_buffer = StringIO()
        read_count = 0
        
        # Paired-end processing with buffered writes
        with dnaio.open(fastq[0], fastq[1], mode="r", fileformat="fastq") as reader, \
             dnaio.open(args.output + "_R1_renamed.fastq.gz", args.output + "_R2_renamed.fastq.gz", 
                      mode="w", fileformat="fastq", compression_level=6) as writer, \
             gzip.open(args.output + "_headers_fw.txt.gz", "wt", encoding="utf-8") as outfile1, \
             gzip.open(args.output + "_headers_rv.txt.gz", "wt", encoding="utf-8") as outfile2:
            
            for record_fw, record_rv in reader:
                header_fw = record_fw.name
                header_rv = record_rv.name
                
                renamed_fw = process_read_header(header_fw, is_forward=True, filename=fastq[0])
                renamed_rv = process_read_header(header_rv, is_forward=False, filename=fastq[1])
                
                # Check if base IDs match
                if renamed_fw != renamed_rv:
                    msg = f"Read IDs were not matching! Please provide matching IDs in the headers. The problematic reads were {header_fw} and {header_rv} in the files {' and '.join(fastq)}."
                    raise ValueError(msg)
                
                fw_buffer.write(f"{renamed_fw}\t{header_fw}\n")
                rv_buffer.write(f"{renamed_rv}\t{header_rv}\n")
                
                record_fw.name = renamed_fw
                record_rv.name = renamed_rv
                
                writer.write(record_fw, record_rv)
                read_count += 1
                
                if read_count >= buffer_size:
                    write_buffer(fw_buffer, outfile1)
                    write_buffer(rv_buffer, outfile2)
                    read_count = 0
            
            # Write any remaining data in buffers
            write_buffer(fw_buffer, outfile1)
            write_buffer(rv_buffer, outfile2)
    else:
        header_buffer = StringIO()
        read_count = 0
        
        # Single-end processing with buffered writes
        with dnaio.open(fastq[0], mode="r", fileformat="fastq") as reader, \
             dnaio.open(args.output + "_renamed.fastq.gz", mode="w", fileformat="fastq", compression_level=6) as writer, \
             gzip.open(args.output + "_headers.txt.gz", "wt", encoding="utf-8") as outfile:
            
            for record in reader:
                header = record.name
                renamed = process_read_header(header, filename=fastq[0])
                
                header_buffer.write(f"{renamed}\t{header}\n")
                record.name = renamed
                writer.write(record)
                read_count += 1
                
                # Write buffer if it reaches the threshold
                if read_count >= buffer_size:
                    write_buffer(header_buffer, outfile)
                    read_count = 0
            
            # Write any remaining data in buffer
            write_buffer(header_buffer, outfile)

if __name__ == "__main__":
    sys.exit(main())

