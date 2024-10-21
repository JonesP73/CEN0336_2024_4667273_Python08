import re  # Importa a biblioteca para usar expressões regulares
import sys  # Importa a biblioteca para lidar com argumentos passados pelo terminal
import os   # Importa a biblioteca para trabalhar com arquivos e caminhos

# Dicionário para armazenar as sequências, onde a chave é o identificador e o valor são as bases A, T, C, G
sequencias = {}

# Verifica se foi passado um arquivo como argumento
if len(sys.argv) < 2:
    print("Nenhum arquivo especificado")
    sys.exit(1)  # Encerra o programa

# Pega o nome do arquivo que foi passado como argumento
arquivo = sys.argv[1]

# Verifica se o arquivo existe
if not os.path.exists(arquivo):
    print("Arquivo não encontrado")
    sys.exit(1)  # Encerra o programa

# Abre o arquivo e lê linha por linha
with open(arquivo) as f:
    for linha in f:
        # Remove espaços em branco no final da linha e converte tudo para maiúsculas
        linha = linha.rstrip().upper()
        
        # Se a linha começar com ">", é um cabeçalho com o identificador
        if linha.startswith(">"):
            id = re.search(r'>(\S+)(\s.+?)', linha)  # Busca o identificador na linha
            identificador = id.group(1)  # Extrai o identificador
            # Cria um dicionário para contar as bases (A, T, C, G) para esse identificador
            sequencias[identificador] = {"A": 0, "T": 0, "C": 0, "G": 0}
        else:
            # Conta quantas vezes cada base aparece e atualiza no dicionário
            sequencias[identificador]["A"] += linha.count("A")
            sequencias[identificador]["T"] += linha.count("T")
            sequencias[identificador]["C"] += linha.count("C")
            sequencias[identificador]["G"] += linha.count("G")

# Exibe os resultados
for id in sequencias:
    print(f"{id}\tA_{sequencias[id]['A']}\tT_{sequencias[id]['T']}\tC_{sequencias[id]['C']}\tG_{sequencias[id]['G']}")
