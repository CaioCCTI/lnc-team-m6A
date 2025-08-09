#!/usr/bin/env python3

import re
from Bio import SeqIO

def find_drach_motifs(fasta_file, output_file=None):
    # Define the DRACH regex pattern
    drach_regex = re.compile(r'[AGU][AG]A[C][ACU]')

    results = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper().replace("T", "U")  # Convert DNA to RNA
        seq_id = record.id

        for match in drach_regex.finditer(seq):
            start = match.start() + 1  # 1-based indexing
            end = match.end()
            motif = match.group()
            results.append((seq_id, start, end, motif))

    # Output to terminal or file
    if output_file:
        with open(output_file, 'w') as out:
            out.write("SeqID\tStart\tEnd\tMotif\n")
            for r in results:
                out.write("\t".join(map(str, r)) + "\n")
    else:
        print("SeqID\tStart\tEnd\tMotif")
        for r in results:
            print("\t".join(map(str, r)))

# Example usage
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Find DRACH motifs in a FASTA file.")
    parser.add_argument("fasta", help="Input FASTA file (DNA or RNA)")
    parser.add_argument("-o", "--output", help="Output file to save results (optional)")
    args = parser.parse_args()

    find_drach_motifs(args.fasta, args.output)

