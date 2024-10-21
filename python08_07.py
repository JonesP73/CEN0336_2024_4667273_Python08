import re
import sys
import os

# Dicionário para armazenar identificadores e suas sequências
sequencias = {}

# Verifica se um arquivo foi passado como argumento
if len(sys.argv) < 2:
    print("Nenhum arquivo fornecido.")
    sys.exit(1)

# Obtém o nome do arquivo do argumento
arquivo = sys.argv[1]

# Verifica se o arquivo existe
if not os.path.exists(arquivo):
    print("Arquivo não encontrado.")
    sys.exit(1)

# Abre o arquivo e processa suas linhas
with open(arquivo) as f:
    identificador = ""
    for linha in f:
        linha = linha.rstrip()  # Remove espaços em branco no final
        if linha.startswith(">"):  # Se for uma linha de cabeçalho
            if identificador:  # Se já existe um identificador, calcula a sequência reversa
                sequencias[identificador]["reverse_seq"] = sequencias[identificador]["seq"][::-1].translate(str.maketrans("ATCG", "TAGC"))
            id = re.search(r'>(\S+)(\s.+?)', linha)  # Extrai o identificador
            identificador = id.group(1)  # Atualiza o identificador
            # Cria um dicionário para armazenar a sequência e os frames
            sequencias[identificador] = {"seq": "", "frame_+1": [], "frame_+2": [], "frame_+3": [], "reverse_seq": "", "frame_-1": [], "frame_-2": [], "frame_-3": []}
        else:  # Se for uma linha de sequência
            sequencias[identificador]["seq"] += linha.upper()  # Adiciona a sequência em maiúsculas

    # Calcula a sequência reversa para o último identificador
    sequencias[identificador]["reverse_seq"] = sequencias[identificador]["seq"][::-1].translate(str.maketrans("ATCG", "TAGC"))

# Abre um novo arquivo para escrever os resultados
with open("Python_08.codons-6frames.nt", "w") as f:
    for id in sequencias:
        # Cria strings de saída
        saida = ""
        saida_rev = ""
        for i in range(3):  # Para cada frame
            # Frame positivo
            frame = f"frame_+{i+1}"
            for match in re.finditer(r"(.{3})", sequencias[id]["seq"][i:]):  # Busca grupos de 3 bases
                sequencias[id][frame].append(match.group(1))  # Adiciona ao dicionário
                if len(sequencias[id][frame][-1]) != 3:  # Remove último elemento se não tiver 3 bases
                    sequencias[id][frame][-1] = sequencias[id][frame][-1][:-1]  # Remove o último caractere

            # Frame reverso
            frame = f"frame_-{i+1}"
            for match in re.finditer(r"(.{3})", sequencias[id]["reverse_seq"][i:]):
                sequencias[id][frame].append(match.group(1))
                if len(sequencias[id][frame][-1]) != 3:
                    sequencias[id][frame][-1] = sequencias[id][frame][-1][:-1]

            # Cria a saída no formato desejado
            saida += f"{id}-frame-{i+1}-codons\n{' '.join(sequencias[id][frame])}\n"
            saida_rev += f"{id}-frame--{i+1}-codons\n{' '.join(sequencias[id][frame])}\n"

        # Escreve os resultados no arquivo
        f.write(saida)
        f.write(saida_rev)

# Tabela de tradução
tabela_traducao = {
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    'AAT': 'N', 'AAC': 'N',
    'GAT': 'D', 'GAC': 'D',
    'TGT': 'C', 'TGC': 'C',
    'CAA': 'Q', 'CAG': 'Q',
    'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    'CAT': 'H', 'CAC': 'H',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
    'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'AAA': 'K', 'AAG': 'K',
    'ATG': 'M',
    'TTT': 'F', 'TTC': 'F',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'TGG': 'W',
    'TAT': 'Y', 'TAC': 'Y',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TAA': '*', 'TGA': '*', 'TAG': '*'
}

# Abre arquivos para salvar as traduções
with open("Python_08.translated.aa", "w") as f:
    with open("Python_08.translated-longest.aa", "w") as longest_f:
        with open("Python_08.orf-longest.nt", "w") as orf_f:
            for id in sequencias:
                maior_peptideo = ""
                maior_peptideo_len = 0
                sequencia_maior_peptideo = ""

                # Processa cada frame positivo
                for i in range(3):
                    frame = f"frame_+{i+1}"
                    protein_frame = f"protein_frame_+{i+1}"
                    sequencias[id][protein_frame] = ""  # Inicializa a proteína do frame

                    for codon in sequencias[id][frame]:
                        if codon not in tabela_traducao:  # Se o códon não estiver na tabela, usa "_"
                            print(f"Aviso: códon {codon} não encontrado, usando _ como substituto.")
                            sequencias[id][protein_frame] += "_"
                        else:
                            sequencias[id][protein_frame] += tabela_traducao[codon]

                    # Busca o maior peptídeo traduzido
                    for peptide in re.finditer(r"(M[A-Z]+?)\*", sequencias[id][protein_frame]):
                        if len(peptide.group(1)) > maior_peptideo_len:
                            maior_peptideo_len = len(peptide.group(1))
                            maior_peptideo = peptide.group(1)
                            sequencia_maior_peptideo = "".join(sequencias[id][frame][peptide.start():peptide.end() + 1])

                # Processa cada frame reverso
                for i in range(3):
                    frame = f"frame_-{i+1}"
                    protein_frame = f"protein_frame_-{i+1}"
                    sequencias[id][protein_frame] = ""

                    for codon in sequencias[id][frame]:
                        if codon not in tabela_traducao:
                            print(f"Aviso: códon {codon} não encontrado, usando _ como substituto.")
                            sequencias[id][protein_frame] += "_"
                        else:
                            sequencias[id][protein_frame] += tabela_traducao[codon]

                    # Busca o maior peptídeo traduzido
                    for peptide in re.finditer(r"(M[A-Z]+?)\*", sequencias[id][protein_frame]):
                        if len(peptide.group(1)) > maior_peptideo_len:
                            maior_peptideo_len = len(peptide.group(1))
                            maior_peptideo = peptide.group(1)
                            sequencia_maior_peptideo = "".join(sequencias[id][frame][peptide.start():peptide.end() + 1])

                # Salva o maior peptídeo traduzido
                longest_f.write(f">{id}\n{maior_peptideo}\n")
                orf_f.write(f">{id}\n{sequencia_maior_peptideo}\n")

                # Salva os resultados em formato FASTA
                f.write("".join([f">{id}-frame-{i}-protein\n{sequencias[id]['protein_frame_+' + str(i)]}\n" for i in range(1, 4)]))
                f.write("".join([f">{id}-frame--{i}-protein\n{sequencias[id]['protein_frame_-' + str(i)]}\n" for i in range(1, 4)]))
