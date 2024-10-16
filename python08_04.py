from Bio import SeqIO
from Bio.Seq import Seq

# Solicita o nome do arquivo FASTA ao usuário
input_file = input("Digite o nome do arquivo FASTA (ex: Python_08.fasta): ")

# Define o nome do arquivo de saída
output_file = "Python_08.codons-6frames.nt"

# Abre o arquivo de saída para escrita
with open(output_file, "w") as outfile:
    # Itera por cada sequência no arquivo FASTA
    for record in SeqIO.parse(input_file, "fasta"):
        # Obtém o ID da sequência
        sequence_id = record.id
        sequence = record.seq
        
        # Processa os três quadros de leitura da sequência original
        for frame in range(3):
            # Divide a sequência em códons no quadro de leitura atual
            codons = [str(sequence[i:i+3]) for i in range(frame, len(sequence), 3) if len(sequence[i:i+3]) == 3]
            
            # Escreve o cabeçalho e os códons no arquivo de saída
            outfile.write(f">{sequence_id}-frame-{frame+1}-codons\n")
            outfile.write(" ".join(codons) + "\n")
        
        # Calcula o complemento reverso da sequência
        reverse_complement = sequence.reverse_complement()
        
        # Processa os três quadros de leitura do complemento reverso
        for frame in range(3):
            # Divide a sequência em códons no quadro de leitura atual
            codons = [str(reverse_complement[i:i+3]) for i in range(frame, len(reverse_complement), 3) if len(reverse_complement[i:i+3]) == 3]
            
            # Escreve o cabeçalho e os códons no arquivo de saída
            outfile.write(f">{sequence_id}-reverse-frame-{frame+1}-codons\n")
            outfile.write(" ".join(codons) + "\n")

print(f"Arquivo '{output_file}' criado com sucesso!")

