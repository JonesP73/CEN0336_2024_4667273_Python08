from Bio import SeqIO
from Bio.Seq import Seq

# Tabela de tradução fornecida
tabela_de_traducao = {
    'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
    'AAT':'N', 'AAC':'N',
    'GAT':'D', 'GAC':'D',
    'TGT':'C', 'TGC':'C',
    'CAA':'Q', 'CAG':'Q',
    'GAA':'E', 'GAG':'E',
    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
    'CAT':'H', 'CAC':'H',
    'ATT':'I', 'ATC':'I', 'ATA':'I',
    'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
    'AAA':'K', 'AAG':'K',
    'ATG':'M',
    'TTT':'F', 'TTC':'F',
    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'TGG':'W',
    'TAT':'Y', 'TAC':'Y',
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
    'TAA':'*', 'TGA':'*', 'TAG':'*'
}

# Solicita o nome do arquivo FASTA ao usuário
input_file = input("Digite o nome do arquivo FASTA (ex: Python_08.fasta): ")

# Define os nomes dos arquivos de saída
output_codons_file = "Python_08.codons-6frames.nt"
output_translation_file = "Python_08.translated.aa"

# Abre os arquivos de saída para escrita
with open(output_codons_file, "w") as codons_outfile, open(output_translation_file, "w") as translation_outfile:
    # Itera por cada sequência no arquivo FASTA
    for record in SeqIO.parse(input_file, "fasta"):
        # Obtém o ID e a sequência da sequência
        sequence_id = record.id
        sequence = record.seq
        
        # Processa os três quadros de leitura da sequência original
        for frame in range(3):
            # Divide a sequência em códons no quadro de leitura atual
            codons = [str(sequence[i:i+3]) for i in range(frame, len(sequence), 3) if len(sequence[i:i+3]) == 3]
            
            # Escreve o cabeçalho e os códons no arquivo de códons
            codons_outfile.write(f">{sequence_id}-frame-{frame+1}-codons\n")
            codons_outfile.write(" ".join(codons) + "\n")
            
            # Traduz os códons em aminoácidos
            amino_acids = [tabela_de_traducao.get(codon, 'X') for codon in codons]
            translation_outfile.write(f">{sequence_id}-frame-{frame+1}-aa\n")
            translation_outfile.write("".join(amino_acids) + "\n")
        
        # Calcula o complemento reverso da sequência
        reverse_complement = sequence.reverse_complement()
        
        # Processa os três quadros de leitura do complemento reverso
        for frame in range(3):
            # Divide a sequência em códons no quadro de leitura atual
            codons = [str(reverse_complement[i:i+3]) for i in range(frame, len(reverse_complement), 3) if len(reverse_complement[i:i+3]) == 3]
            
            # Escreve o cabeçalho e os códons no arquivo de códons
            codons_outfile.write(f">{sequence_id}-reverse-frame-{frame+1}-codons\n")
            codons_outfile.write(" ".join(codons) + "\n")
            
            # Traduz os códons em aminoácidos
            amino_acids = [tabela_de_traducao.get(codon, 'X') for codon in codons]
            translation_outfile.write(f">{sequence_id}-reverse-frame-{frame+1}-aa\n")
            translation_outfile.write("".join(amino_acids) + "\n")

print(f"Arquivos '{output_codons_file}' e '{output_translation_file}' criados com sucesso!")
