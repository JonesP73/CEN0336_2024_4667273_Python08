import re
import sys
import os

# Dicionário para armazenar as sequências e informações associadas a cada identificador
seqs = {}

# Verifica se um arquivo foi passado como argumento
if len(sys.argv) < 2:
    print("Nenhum arquivo especificado")
    sys.exit(1)

arquivo = sys.argv[1]

# Verifica se o arquivo existe
if not os.path.exists(arquivo):
    print("Arquivo não encontrado")
    sys.exit(1)

# Lê o arquivo e processa as sequências
with open(arquivo) as f:
    identificador = ""
    for linha in f:
        linha = linha.rstrip()  # Remove espaços e quebras de linha
        if linha.startswith(">"):  # Verifica se é um cabeçalho
            if identificador:
                # Cria a sequência complementar reversa
                seqs[identificador]["rev_comp_seq"] = seqs[identificador]["seq"][::-1].translate(str.maketrans("ATCG", "TAGC"))
            
            # Extrai o identificador
            id_match = re.search(r'>(\S+)', linha)
            identificador = id_match.group(1)
            # Cria um dicionário para cada identificador
            seqs[identificador] = {
                "seq": "",
                "frame_+1": [], "frame_+2": [], "frame_+3": [],
                "rev_comp_seq": "",
                "frame_-1": [], "frame_-2": [], "frame_-3": []
            }
        else:
            # Adiciona a sequência ao dicionário e coloca em maiúsculas
            seqs[identificador]["seq"] += linha.upper()

    # Cria a sequência complementar reversa para o último identificador
    seqs[identificador]["rev_comp_seq"] = seqs[identificador]["seq"][::-1].translate(str.maketrans("ATCG", "TAGC"))

# Abre o arquivo para salvar os códons em 6 frames
with open("Python_08.codons-6frames.nt", "w") as f:
    for id in seqs:
        output = ""
        output_rev = ""
        for i in range(3):
            # Frame positivo
            frame = f"frame_+{i+1}"
            for match in re.finditer(r"(.{3})", seqs[id]["seq"][i:]):
                seqs[id][frame].append(match.group(1))

            # Frame negativo
            frame_rev = f"frame_-{i+1}"
            for match in re.finditer(r"(.{3})", seqs[id]["rev_comp_seq"][i:]):
                seqs[id][frame_rev].append(match.group(1))

        # Escreve os resultados no arquivo
        for i in range(3):
            output += f"{id}-frame-{i+1}-codons\n{' '.join(seqs[id][f'frame_+{i+1}'])}\n"
            output_rev += f"{id}-frame--{i+1}-codons\n{' '.join(seqs[id][f'frame_-{i+1}'])}\n"
        f.write(output)
        f.write(output_rev)

# Tabela de tradução de códons para aminoácidos
tabela_traducao = {
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

# Abre o arquivo para salvar as traduções em proteína
with open("Python_08.translated.aa", "w") as f:
    for id in seqs:
        for i in range(3):
            frame = f"frame_+{i+1}"
            protein_frame = f"protein_frame_+{i+1}"
            seqs[id][protein_frame] = "".join(
                tabela_traducao.get(codon, "_") for codon in seqs[id][frame]
            )

            frame_rev = f"frame_-{i+1}"
            protein_frame_rev = f"protein_frame_-{i+1}"
            seqs[id][protein_frame_rev] = "".join(
                tabela_traducao.get(codon, "_") for codon in seqs[id][frame_rev]
            )

        # Escreve as traduções no arquivo
        for i in range(1, 4):
            f.write(f">{id}-frame-{i}-protein\n{seqs[id][f'protein_frame_+{i}']}\n")
            f.write(f">{id}-frame--{i}-protein\n{seqs[id][f'protein_frame_-{i}']}\n")
