#!/usr/bin/env python3

# Função para ler um arquivo FASTA e retornar um dicionário com sequências
def read_fasta(file_path):
    """Lê um arquivo FASTA e retorna um dicionário com sequências."""
    sequences = {}
    seq_id = ""
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                seq_id = line[1:]  # Remove o caractere '>'
                sequences[seq_id] = ""
            else:
                sequences[seq_id] += line  # Acumula a sequência

    return sequences

# Função para gerar códons a partir de uma sequência
def get_codons(sequence, frame):
    """Divide a sequência em códons (trincas de nucleotídeos) no quadro especificado."""
    sequence = sequence[frame:]  # Ajusta a sequência para o quadro de leitura
    return [sequence[i:i + 3] for i in range(0, len(sequence) - 2, 3)]  # Divide em códons de 3 em 3

# Função para obter o complemento reverso de uma sequência
def reverse_complement(sequence):
    """Retorna o complemento reverso da sequência."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}  # Dicionário de complementos
    rev_comp = ''.join(complement[base] for base in reversed(sequence))  # Gera o complemento reverso
    return rev_comp

# Solicita o caminho do arquivo FASTA ao usuário
fasta_file = input("Digite o caminho completo do arquivo Python_08.codons-3frames.nt: ")
output_file = 'Python_08.codons-6frames.nt'

try:
    with open(output_file, 'w') as out_file:
        # Lê as sequências do arquivo FASTA
        sequences = read_fasta(fasta_file)
        
        for seq_id, sequence in sequences.items():
            # Obter códons para os três quadros de leitura da sequência original
            for frame in range(3):  # Para os três quadros de leitura
                codons = get_codons(sequence, frame)
                out_file.write(f"{seq_id}-frame-{frame + 1}-codons\n")
                out_file.write(" ".join(codons) + "\n")

            # Obter o complemento reverso da sequência
            rev_comp = reverse_complement(sequence)
            # Obter códons para os três quadros de leitura do complemento reverso
            for frame in range(3):
                codons = get_codons(rev_comp, frame)
                out_file.write(f"{seq_id}-rev-comp-frame-{frame + 1}-codons\n")
                out_file.write(" ".join(codons) + "\n")

    print(f"Códons foram escritos em {output_file}.")
except FileNotFoundError:
    print(f"Erro: O arquivo {fasta_file} não foi encontrado. Verifique o caminho e tente novamente.")
except Exception as e:
    print(f"Um erro ocorreu: {e}")
