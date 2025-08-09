#!/usr/bin/env python3

import sys

def recorta_fasta_dna_para_rna(arquivo_entrada, arquivo_saida, inicio, fim):
    with open(arquivo_entrada, 'r') as fin, open(arquivo_saida, 'w') as fout:
        seq_id = ''
        seq = ''
        for line in fin:
            line = line.strip()
            if line.startswith('>'):
                if seq_id:
                    # Aplica corte + substituição T > U
                    recortada = cortar_e_converter(seq, inicio, fim)
                    fout.write(f"{seq_id}\n{recortada}\n")
                seq_id = line
                seq = ''
            else:
                seq += line
        if seq_id:
            recortada = cortar_e_converter(seq, inicio, fim)
            fout.write(f"{seq_id}\n{recortada}\n")

def cortar_e_converter(seq, inicio, fim):
    # Corrige índices para base 0
    if inicio < fim:
        trecho = seq[inicio-1:fim]
    else:
        trecho = seq[fim-1:inicio][::-1]  # reverso se fim < início
    # Substitui T/t por U/u
    return trecho.replace('T', 'U').replace('t', 'u')

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Uso: python recorta_fasta.py entrada.fasta saida.fasta inicio fim")
        sys.exit(1)

    entrada = sys.argv[1]
    saida = sys.argv[2]
    inicio = int(sys.argv[3])
    fim = int(sys.argv[4])

    recorta_fasta_dna_para_rna(entrada, saida, inicio, fim)

