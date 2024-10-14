#!/usr/bin/env python3

# Função para gerar códons a partir de uma sequência
def get_codons(sequence, frame):
    # Ajusta a sequência para o quadro de leitura
    sequence = sequence[frame:]
    # Divide a sequência em códons de 3 em 3 e retorna como uma lista
    return [sequence[i:i + 3] for i in range(0, len(sequence) - 2, 3)]

# Solicita o caminho do arquivo de entrada ao usuário
codons_file = input("Digite o caminho completo do arquivo Python_08.codons-frame-1.nt: ")
output_file = 'Python_08.codons-3frames.nt'

try:
    with open(codons_file, 'r') as in_file, open(output_file, 'w') as out_file:
        # Lê as sequências do arquivo de códons
        for line in in_file:
            line = line.strip()
            if line.endswith('-codons'):
                record_id = line[:-8]  # Remove '-codons' para obter o ID
                # Lê a linha seguinte que contém os códons
                codons_line = next(in_file).strip()
                codons = codons_line.split()  # Divide os códons em uma lista

                # Escreve os resultados no arquivo de saída para os três quadros de leitura
                full_sequence = "".join(codons)  # Junta os códons em uma sequência
                for frame in range(3):  # Para os três quadros de leitura
                    frame_codons = get_codons(full_sequence, frame)
                    out_file.write(f"{record_id}-frame-{frame + 1}-codons\n")
                    out_file.write(" ".join(frame_codons) + "\n")

    print(f"Códons foram escritos em {output_file}.")
except FileNotFoundError:
    print(f"Erro: O arquivo {codons_file} não foi encontrado. Verifique o caminho e tente novamente.")
except Exception as e:
    print(f"Um erro ocorreu: {e}")
