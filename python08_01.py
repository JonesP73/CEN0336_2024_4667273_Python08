def read_fasta(file_path):
    # Dicionário para armazenar as sequências
    seqs = {}
    gene_name = ""

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Linha com nome da sequência
                gene_name = line[1:]  # Remove o caractere '>'
                seqs[gene_name] = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
            else:
                # Linha com a sequência (acumulando contagem de nucleotídeos)
                for nucleotide in line:
                    if nucleotide in seqs[gene_name]:
                        seqs[gene_name][nucleotide] += 1

    return seqs

def print_nucleotide_composition(seqs):
    for gene_name, counts in seqs.items():
        print(f"{gene_name}\t{counts['A']}\t{counts['T']}\t{counts['G']}\t{counts['C']}")  # Fixed syntax error: Added closing bracket

# Solicita o caminho do arquivo fornecido pelo usuário
file_path = input("Digite o caminho do arquivo multi-FASTA: ")

# Lendo e processando o arquivo
sequences = read_fasta(file_path)

# Imprimindo a composição de nucleotídeos
print_nucleotide_composition(sequences)
