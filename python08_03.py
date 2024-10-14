from Bio import SeqIO

# Função para gerar códons a partir de uma sequência
def get_codons(sequence, frame):
    codons = []
    # Ajusta a sequência para o quadro de leitura
    sequence = sequence[frame:]
    # Divide a sequência em códons de 3 em 3
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i + 3]
        codons.append(codon)
    return codons

# Solicita o caminho do arquivo FASTA ao usuário
fasta_file = input("Digite o caminho completo do arquivo FASTA (ex: /caminho/para/Python_08.fasta): ")
output_file = 'Python_08.codons-3frames.nt'

try:
    with open(output_file, 'w') as out_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            for frame in range(3):  # Para os três quadros de leitura
                codons = get_codons(str(record.seq), frame)
                # Escreve os resultados no arquivo
                out_file.write(f"{record.id}-frame-{frame+1}-codons\n")
                out_file.write(" ".join(codons) + "\n")

    print(f"Códons foram escritos em {output_file}.")
except FileNotFoundError:
    print(f"Erro: O arquivo {fasta_file} não foi encontrado. Verifique o caminho e tente novamente.")
except Exception as e:
    print(f"Um erro ocorreu: {e}")
