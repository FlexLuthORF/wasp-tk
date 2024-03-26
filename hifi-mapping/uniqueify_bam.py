#!/usr/bin/env python

import pysam
import subprocess
import sys

def uniqueify_bam_read_names(input_bam, output_bam):
    # Open the input BAM file for reading
    with pysam.AlignmentFile(input_bam, "rb") as infile:
        # Open a new BAM file for writing with the same header as the input
        with pysam.AlignmentFile(output_bam, "wb", header=infile.header) as outfile:
            # Iterate over each read in the input BAM file
            for i, read in enumerate(infile):
                # Modify the read name to ensure uniqueness
                new_read_name = f"{read.query_name}_{i}"
                read.query_name = new_read_name
                
                # Write the modified read to the output BAM file
                outfile.write(read)
    
    # Index the output BAM file using samtools
    subprocess.run(["samtools", "index", output_bam])

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python uniqueify_bam.py <input.bam> <output.bam>", file=sys.stderr)
        sys.exit(1)

    input_bam_path = sys.argv[1]
    output_bam_path = sys.argv[2]

    uniqueify_bam_read_names(input_bam_path, output_bam_path)
