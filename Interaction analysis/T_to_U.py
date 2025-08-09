#!/usr/bin/env python3

import sys

def dna_to_rna_fasta(input_file, output_file):
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                fout.write(line)
            else:
                rna_line = line.strip().replace('T', 'U').replace('t', 'u')
                fout.write(rna_line + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Uso: python dna_to_rna_fasta.py entrada.fasta saida.fasta")
        sys.exit(1)

    entrada = sys.argv[1]
    saida = sys.argv[2]
    dna_to_rna_fasta(entrada, saida)

