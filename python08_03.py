import re  # Importa o módulo para trabalhar com expressões regulares
import sys  # Importa o módulo para lidar com argumentos da linha de comando
import os   # Importa o módulo para trabalhar com arquivos

# Dicionário para armazenar identificadores e suas sequências
sequencias = {}

# Verifica se um arquivo foi passado como argumento
if len(sys.argv) < 2:
    print("Nenhum arquivo especificado")
    sys.exit(1)  # Encerra o programa se não houver arquivo

# Pega o nome do arquivo que foi passado como argumento
arquivo = sys.argv[1]

# Verifica se o arquivo existe
if not os.path.exists(arquivo):
    print("Arquivo não encontrado")
    sys.exit(1)  # Encerra o programa se o arquivo não existir

# Abre o arquivo e lê cada linha
with open(arquivo) as f:
    for linha in f:
        linha = linha.rstrip()  # Remove espaços e quebras de linha no final

        # Verifica se a linha é um cabeçalho (começa com ">")
        if linha.startswith(">"):
            id = re.search(r'>(\S+)(\s.+?)', linha)  # Busca o identificador
            identificador = id.group(1)  # Pega o identificador encontrado
            # Cria um dicionário para armazenar a sequência e os três frames
            sequencias[identificador] = {"seq": "", "frame_+1": [], "frame_+2": [], "frame_+3": []}
        else:
            # Adiciona a sequência convertida para maiúsculas ao identificador atual
            sequencias[identificador]["seq"] += linha.upper()

# Abre o arquivo de saída para escrever os resultados
with open("Python_08.codons-3frames.nt", "w") as f:
    for id in sequencias:  # Para cada identificador no dicionário
        # Itera sobre os três frames (0, 1 e 2)
        for i in range(3):
            frame = "frame_+" + str(i + 1)
            # Busca grupos de 3 bases na sequência a partir do índice do frame
            for match in re.finditer(r"(.{3})", sequencias[id]["seq"][i:]):
                sequencias[id][frame].append(match.group(1))  # Adiciona o grupo de 3 ao frame

            # Remove o último grupo se ele não tiver 3 bases
            if len(sequencias[id][frame][-1]) != 3:
                sequencias[id][frame].pop()

            # Escreve os resultados formatados no arquivo de saída
            f.write(f"{id}-frame-{i+1}-codons\n{' '.join(sequencias[id][frame])}\n")
