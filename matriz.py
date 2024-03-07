from functions import getKmersSeq, assign_family, assign_family_or_op
import pandas as pd
import json


sequencia = input("Informe a sequência a ser classificada\n").upper()
k_cl = getKmersSeq(sequencia, 8)



# lendo k_mers do database PR2
with open('teste.json') as f:
    k_db = json.load(f)


# gerando classificao 
classificacao = assign_family(k_cl, k_db)
classificacao_or = assign_family_or_op(classificacao, k_cl)

# escrevendo matriz familia op or
df = (pd.DataFrame(classificacao_or)).transpose().sort_index(axis=1)
df.to_excel('teste.xlsx', index=True)