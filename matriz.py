from collections import defaultdict
from functions import getKmersSeq,getKmers, assign_family, assign_family_or_op
import pandas as pd
import json

# leitura do database
arq = open("pr2_sample.fasta", "r")
linhas = arq.readlines()
organismos = defaultdict(int)
count = defaultdict(int)
all_families = set()
for i in range(0, len(linhas), 2):
  all_families.add(linhas[i].split(";")[6])
  count[linhas[i].split(";")[6]] +=1
  organismos[linhas[i].split(";")[6] + " " + str(count[linhas[i].split(";")[6]])] = linhas[i+1]

sequencia = input("Informe a sequência a ser classificada\n").upper()

k_cl = getKmersSeq(sequencia, 8)
k_db = getKmers(organismos, 8)

# # lendo k_mers do database PR2 (caso haja um json com o database pré-processado)
# with open('teste.json') as f:
#     k_db = json.load(f)


# gerando classificao 
classificacao = assign_family(k_cl, k_db)
classificacao_or = assign_family_or_op(classificacao, k_cl, all_families)

# escrevendo matriz familia op or
df = (pd.DataFrame(classificacao_or)).transpose().sort_index(axis=1)
df.to_excel('teste.xlsx', index=True)