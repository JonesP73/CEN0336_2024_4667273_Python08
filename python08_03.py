from Bio import SeqIO

# Solicita o nome do arquivo FASTA ao usuário
input_file = input("Digite o nome do arquivo FASTA (ex: Python_08.fasta): ")

# Define o nome do arquivo de saída
output_file = "Python_08.codons-3frames.nt"

# Abre o arquivo de saída para escrita
with open(output_file, "w") as outfile:
    # Itera por cada sequência no arquivo FASTA
    for record in SeqIO.parse(input_file, "fasta"):
        # Obtém o ID da sequência
        sequence_id = record.id
        
        # Itera pelos três primeiros quadros de leitura
        for frame in range(3):
            # Divide a sequência em códons no quadro de leitura atual
            codons = [str(record.seq[i:i+3]) for i in range(frame, len(record.seq), 3) if len(record.seq[i:i+3]) == 3]
            
            # Escreve o cabeçalho e os códons no arquivo de saída
            outfile.write(f">{sequence_id}-frame-{frame+1}-codons\n")
            outfile.write(" ".join(codons) + "\n")

print(f"Arquivo '{output_file}' criado com sucesso!")
