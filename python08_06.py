import re
import sys
import os

# Dicionário para armazenar os identificadores como chave e as sequências como valor
seqs = {}

# Verifica se um arquivo foi fornecido como argumento
if len(sys.argv) < 2:
    print("Nenhum arquivo especificado.")
    sys.exit(1)

file = sys.argv[1]

# Verifica se o arquivo existe
if not os.path.exists(file):
    print("Arquivo não encontrado.")
    sys.exit(1)

# Abrindo o arquivo e processando linha por linha
with open(file) as f:
    identificador = ""
    for linha in f:
        linha = linha.rstrip()  # Remove espaços e quebras de linha extras
        if linha.startswith(">"):  # Verifica se é um cabeçalho
            # Se já existe um identificador, cria a sequência reversa complementar
            if identificador:
                seqs[identificador]["reverse_seq"] = seqs[identificador]["seq"][::-1].translate(str.maketrans("ATCG", "TAGC"))
            
            # Extrai o identificador do cabeçalho
            id = re.search(r'>(\S+)(\s.+?)', linha)
            identificador = id.group(1)
            # Cria um dicionário para armazenar as informações da sequência
            seqs[identificador] = {
                "seq": "", 
                "frame_+1": [], "frame_+2": [], "frame_+3": [], 
                "reverse_seq": "", 
                "frame_-1": [], "frame_-2": [], "frame_-3": []
            }
        else:
            # Adiciona a sequência de bases ao dicionário
            seqs[identificador]["seq"] += linha.upper()

    # Cria a sequência reversa complementar para o último identificador
    seqs[identificador]["reverse_seq"] = seqs[identificador]["seq"][::-1].translate(str.maketrans("ATCG", "TAGC"))

# Escreve as sequências em 6 frames em um arquivo
with open("Python_08.codons-6frames.nt", "w") as f:
    for id in seqs:
        saida = ""
        saida_reversa = ""
        for i in range(3):  # Percorre os 3 frames de leitura
            frame = f"frame_+{i+1}"
            for match in re.finditer(r"(.{3})", seqs[id]["seq"][i:]):
                seqs[id][frame].append(match.group(1))

            if len(seqs[id][frame][-1]) != 3:
                seqs[id][frame].pop()

            # Frames reversos
            frame_reverso = f"frame_-{i+1}"
            for match in re.finditer(r"(.{3})", seqs[id]["reverse_seq"][i:]):
                seqs[id][frame_reverso].append(match.group(1))

            if len(seqs[id][frame_reverso][-1]) != 3:
                seqs[id][frame_reverso].pop()

            # Formata as saídas para salvar
            saida += f">{id}-frame-{i+1}-codons\n{' '.join(seqs[id][frame])}\n"
            saida_reversa += f">{id}-frame--{i+1}-codons\n{' '.join(seqs[id][frame_reverso])}\n"

        f.write(saida)
        f.write(saida_reversa)

# Tabela de tradução para aminoácidos
translation_table = {
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

# Traduz as sequências de nucleotídeos para proteínas e salva os resultados
with open("Python_08.translated.aa", "w") as f:
    with open("Python_08.translated-longest.aa", "w") as longest_f:
        for id in seqs:
            maior_peptideo = ""
            tamanho_maior_peptideo = 0
            for i in range(3):
                frame = f"frame_+{i+1}"
                protein_frame = f"protein_frame_+{i+1}"
                seqs[id][protein_frame] = ""

                for codon in seqs[id][frame]:
                    if codon not in translation_table:
                        seqs[id][protein_frame] += "_"
                    else:
                        seqs[id][protein_frame] += translation_table[codon]

                peptideos = re.findall(r"(M[A-Z]+?)\*", seqs[id][protein_frame])
                maior_peptideo_frame = max(peptideos, key=len) if peptideos else ""

                if len(maior_peptideo_frame) > tamanho_maior_peptideo:
                    tamanho_maior_peptideo = len(maior_peptideo_frame)
                    maior_peptideo = maior_peptideo_frame

                # Faz a mesma coisa para os frames negativos
                frame_reverso = f"frame_-{i+1}"
                protein_frame_reverso = f"protein_frame_-{i+1}"
                seqs[id][protein_frame_reverso] = ""

                for codon in seqs[id][frame_reverso]:
                    if codon not in translation_table:
                        seqs[id][protein_frame_reverso] += "_"
                    else:
                        seqs[id][protein_frame_reverso] += translation_table[codon]

                peptideos_reversos = re.findall(r"(M[A-Z]+?)\*", seqs[id][protein_frame_reverso])
                maior_peptideo_frame_reverso = max(peptideos_reversos, key=len) if peptideos_reversos else ""

                if len(maior_peptideo_frame_reverso) > tamanho_maior_peptideo:
                    tamanho_maior_peptideo = len(maior_peptideo_frame_reverso)
                    maior_peptideo = maior_peptideo_frame_reverso

            longest_f.write(f">{id}\n{maior_peptideo}\n")
            f.write("".join([f">{id}-frame-{i}-protein\n{seqs[id][f'protein_frame_+{i}']}\n" for i in range(1, 4)]))
            f.write("".join([f">{id}-frame--{i}-protein\n{seqs[id][f'protein_frame_-{i}']}\n" for i in range(1, 4)]))
