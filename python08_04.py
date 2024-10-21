import re
import sys
import os

seqs = {} # dicionario para armazenar os identificadores como chave e as bases como valor

if len(sys.argv) < 2:
    print("No file")
    sys.exit(1)

file = sys.argv[1]

if not os.path.exists(file):
    print("File not found")
    sys.exit(1)

with open(file) as f:
    identificador = ""
    for l in f:
        l = l.rstrip()
        if l.startswith(">"):  # testa se é cabeçalho
            if identificador:
                seqs[identificador]["reverse_seq"] = seqs[identificador]["seq"][::-1].translate(str.maketrans("ATCG", "TAGC"))
            id = re.search(r'>(\S+)(\s.+?)', l) # busca o identificador
            identificador = id.group(1)
            seqs[identificador] = {"seq": "", "frame_+1":[], "frame_+2":[], "frame_+3":[], "reverse_seq": "", "frame_-1":[], "frame_-2":[], "frame_-3":[]} # cria um dicionario para o identificador
        else:
            seqs[identificador]["seq"] += l.upper() # adiciona a sequencia ao dicionario
    seqs[identificador]["reverse_seq"] = seqs[identificador]["seq"][::-1].translate(str.maketrans("ATCG", "TAGC"))
    
with open("Python_08.codons-6frames.nt", "w") as f:

    for id in seqs: # para cada identificador
        out = ""
        out_rev = ""

        for i in range(3): # para cada frame, no caso só o +1
        
            frame = "frame_+" + str(i+1) 

            for match in re.finditer(r"(.{3})", seqs[id]["seq"][i:]): # busca todas as bases em grupos de 3
                seqs[id][frame].append(match.group(1)) # adiciona ao dicionario
                
            if len(seqs[id][frame][-1]) != 3: # o último elemento pode ter menos de 3 bases, então removemos se for o caso
                seqs[id][frame][-1].pop() # remove o ultimo elemento
            out += f"{id}-frame-{i+1}-codons\n{' '.join(seqs[id][frame])}\n"
            
            ##Reverse

            frame = "frame_-" + str(i+1) 

            for match in re.finditer(r"(.{3})", seqs[id]["reverse_seq"][i:]): # busca todas as bases em grupos de 3
                seqs[id][frame].append(match.group(1)) # adiciona ao dicionario
                
            if len(seqs[id][frame][-1]) != 3: # o último elemento pode ter menos de 3 bases, então removemos se for o caso
                seqs[id][frame][-1].pop() # remove o ultimo elemento
            
            out_rev += f"{id}-frame--{i+1}-codons\n{' '.join(seqs[id][frame])}\n"
        
        f.write(out)
        f.write(out_rev)
import re
import sys
import os

# Dicionário para armazenar identificadores como chave e as bases (sequências e frames) como valores
seqs = {}

# Verifica se o arquivo foi fornecido como argumento
if len(sys.argv) < 2:
    print("Arquivo não especificado")
    sys.exit(1)

# Pega o nome do arquivo
file = sys.argv[1]

# Verifica se o arquivo existe
if not os.path.exists(file):
    print("Arquivo não encontrado")
    sys.exit(1)

# Abre e lê o arquivo
with open(file) as f:
    identificador = ""
    for linha in f:
        linha = linha.strip()  # Remove espaços e quebras de linha
        if linha.startswith(">"):  # Verifica se a linha é um cabeçalho
            if identificador:
                # Calcula a sequência reversa complementar quando encontramos um novo identificador
                seqs[identificador]["reverse_seq"] = seqs[identificador]["seq"][::-1].translate(str.maketrans("ATCG", "TAGC"))
            
            # Extrai o identificador
            id_match = re.search(r'>(\S+)(\s.+?)', linha)
            identificador = id_match.group(1)

            # Cria um dicionário para armazenar a sequência e os frames para este identificador
            seqs[identificador] = {
                "seq": "",
                "frame_+1": [], "frame_+2": [], "frame_+3": [],
                "reverse_seq": "",
                "frame_-1": [], "frame_-2": [], "frame_-3": []
            }
        else:
            # Adiciona a sequência de bases ao dicionário
            seqs[identificador]["seq"] += linha.upper()

    # Calcula a sequência reversa complementar para o último identificador
    seqs[identificador]["reverse_seq"] = seqs[identificador]["seq"][::-1].translate(str.maketrans("ATCG", "TAGC"))

# Abre o arquivo para salvar os resultados dos frames
with open("Python_08.codons-6frames.nt", "w") as f:
    for id in seqs:  # Para cada identificador
        out = ""
        out_rev = ""
        
        # Processa os frames +1, +2, +3
        for i in range(3):
            frame = "frame_+" + str(i+1)
            # Divide a sequência em grupos de 3 bases
            for match in re.finditer(r"(.{3})", seqs[id]["seq"][i:]):
                seqs[id][frame].append(match.group(1))
                
            # Remove o último elemento se tiver menos de 3 bases
            if len(seqs[id][frame][-1]) != 3:
                seqs[id][frame].pop()
                
            out += f"{id}-frame-{i+1}-codons\n{' '.join(seqs[id][frame])}\n"

        # Processa os frames -1, -2, -3 (reverso)
        for i in range(3):
            frame = "frame_-" + str(i+1)
            for match in re.finditer(r"(.{3})", seqs[id]["reverse_seq"][i:]):
                seqs[id][frame].append(match.group(1))
                
            if len(seqs[id][frame][-1]) != 3:
                seqs[id][frame].pop()
                
            out_rev += f"{id}-frame--{i+1}-codons\n{' '.join(seqs[id][frame])}\n"

        # Escreve os resultados no arquivo
        f.write(out)
        f.write(out_rev)
