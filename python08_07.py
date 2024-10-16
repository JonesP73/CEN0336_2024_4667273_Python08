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

def traduzir_codons(codons):
    return [tabela_de_traducao.get(codon, 'X') for codon in codons]

def encontrar_peptideo_mais_longo(amino_acids, codons):
    longest_peptide = ""
    longest_codons = []
    current_peptide = []
    current_codons = []

    for i, aa in enumerate(amino_acids):
        if aa == 'M':  # Início de um novo peptídeo
            current_peptide = ['M']
            current_codons = [codons[i]]
        elif aa == '*':  # Fim do peptídeo
            if len(current_peptide) > len(longest_peptide):
                longest_peptide = "".join(current_peptide)
                longest_codons = current_codons[:]
            current_peptide = []
            current_codons = []
        elif current_peptide:
            current_peptide.append(aa)
            current_codons.append(codons[i])

    # Verifica se o peptídeo final é o mais longo
    if len(current_peptide) > len(longest_peptide):
        longest_peptide = "".join(current_peptide)
        longest_codons = current_codons

    return longest_peptide, longest_codons

# Solicita o nome do arquivo FASTA ao usuário
input_file = input("Digite o nome do arquivo FASTA (ex: Python_08.fasta): ")

# Define os nomes dos arquivos de saída
output_codons_file = "Python_08.codons-6frames.nt"
output_translation_file = "Python_08.translated.aa"
output_longest_peptide_file = "Python_08.translated-longest.aa"
output_orf_file = "Python_08.orf-longest.nt"

# Abre os arquivos de saída para escrita
with open(output_codons_file, "w") as codons_outfile, \
     open(output_translation_file, "w") as translation_outfile, \
     open(output_longest_peptide_file, "w") as longest_peptide_outfile, \
     open(output_orf_file, "w") as orf_outfile:
    
    # Itera por cada sequência no arquivo FASTA
    for record in SeqIO.parse(input_file, "fasta"):
        sequence_id = record.id
        sequence = record.seq
        all_peptides = []
        all_codons_sets = []

        # Processa os três quadros de leitura da sequência original
        for frame in range(3):
            codons = [str(sequence[i:i+3]) for i in range(frame, len(sequence), 3) if len(sequence[i:i+3]) == 3]
            amino_acids = traduzir_codons(codons)
            
            # Escreve no arquivo de códons
            codons_outfile.write(f">{sequence_id}-frame-{frame+1}-codons\n")
            codons_outfile.write(" ".join(codons) + "\n")
            
            # Escreve no arquivo de tradução
            translation_outfile.write(f">{sequence_id}-frame-{frame+1}-aa\n")
            translation_outfile.write("".join(amino_acids) + "\n")
            
            # Encontra o peptídeo mais longo e seus códons
            peptide, codon_set = encontrar_peptideo_mais_longo(amino_acids, codons)
            all_peptides.append(peptide)
            all_codons_sets.append(codon_set)

        # Processa os três quadros de leitura do complemento reverso
        reverse_complement = sequence.reverse_complement()
        for frame in range(3):
            codons = [str(reverse_complement[i:i+3]) for i in range(frame, len(reverse_complement), 3) if len(reverse_complement[i:i+3]) == 3]
            amino_acids = traduzir_codons(codons)
            
            # Escreve no arquivo de códons
            codons_outfile.write(f">{sequence_id}-reverse-frame-{frame+1}-codons\n")
            codons_outfile.write(" ".join(codons) + "\n")
            
            # Escreve no arquivo de tradução
            translation_outfile.write(f">{sequence_id}-reverse-frame-{frame+1}-aa\n")
            translation_outfile.write("".join(amino_acids) + "\n")
            
            # Encontra o peptídeo mais longo e seus códons
            peptide, codon_set = encontrar_peptideo_mais_longo(amino_acids, codons)
            all_peptides.append(peptide)
            all_codons_sets.append(codon_set)

        # Identifica o peptídeo mais longo e o conjunto de códons correspondente
        longest_peptide = ""
        longest_codons = []
        for i, peptide in enumerate(all_peptides):
            if len(peptide) > len(longest_peptide):
                longest_peptide = peptide
                longest_codons = all_codons_sets[i]

        # Escreve o peptídeo mais longo no arquivo de saída
        longest_peptide_outfile.write(f">{sequence_id}-longest-peptide\n")
        longest_peptide_outfile.write(longest_peptide + "\n")

        # Escreve o conjunto de códons que produziu o peptídeo mais longo
        orf_outfile.write(f">{sequence_id}-longest-orfs\n")
        orf_outfile.write(" ".join(longest_codons) + "\n")

print(f"Arquivos '{output_codons_file}', '{output_translation_file}', '{output_longest_peptide_file}', e '{output_orf_file}' criados com sucesso!")
