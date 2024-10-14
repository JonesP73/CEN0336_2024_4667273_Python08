#!/usr/bin/env python3

def read_fasta(file_path):
    """Lê o arquivo multi-FASTA e retorna um dicionário com sequências."""
    seqs = {}
    gene_name = ""

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Nome da sequência
                gene_name = line[1:]  # Remove o caractere '>'
                seqs[gene_name] = ""
            else:
                # Acumula as sequências
                seqs[gene_name] += line

    return seqs

def get_codons(sequence):
    """Divide a sequência em códons (trincas de nucleotídeos) no primeiro quadro de leitura."""
    return [sequence[i:i + 3] for i in range(0, len(sequence) - len(sequence) % 3, 3)]

def print_codons(seqs):
    """Imprime os códons para cada sequência."""
    for gene_name, sequence in seqs.items():
        codons = get_codons(sequence)
        codons_str = ' '.join(codons)
        print(f"{gene_name}-frame-1-codons")
        print(codons_str)

# Solicita o caminho do arquivo fornecido pelo usuário
file_path = input("Digite o caminho do arquivo multi-FASTA: ")

# Lendo e processando o arquivo
sequences = read_fasta(file_path)

# Imprimindo os códons
print_codons(sequences)

