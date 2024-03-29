from collections import defaultdict
from functions import getKmersSeq,getKmers, assign_family, assign_family_or_op
import pandas as pd
import matplotlib.pyplot as plt

# leitura do database
arq = open("database/pr2_version_5.0.0_SSU_dada2.fasta", "r")
linhas = arq.readlines()
organismos = defaultdict(int)
count = defaultdict(int)
all_families = set()

dict = {
          "Domain":[],
          "Supergroup": [],
          "Division": [],
          "Subdivision": [],
          "Class": [],
          "Order": [],
          "Family": [],
          "Genus": [],
          "Species": [],
          "Sequence": []
        }
      


for i in range(0, len(linhas), 2):
  dict.get("Domain").append(linhas[i].split(";")[0])
  dict.get("Supergroup").append(linhas[i].split(";")[1])
  dict.get("Division").append(linhas[i].split(";")[2])
  dict.get("Subdivision").append(linhas[i].split(";")[3])
  dict.get("Class").append(linhas[i].split(";")[4])
  dict.get("Order").append(linhas[i].split(";")[5])
  dict.get("Family").append(linhas[i].split(";")[6])
  dict.get("Genus").append(linhas[i].split(";")[7])
  dict.get("Species").append(linhas[i].split(";")[8])
  dict.get("Sequence").append(linhas[i+1])

df = pd.DataFrame(dict)

  # all_families.add(linhas[i].split(";")[6])
  # count[linhas[i].split(";")[6]] +=1
  # organismos[linhas[i].split(";")[6] + " " + str(count[linhas[i].split(";")[6]])] = linhas[i+1]
df_heterobranchia = df.loc[df['Family']=="Heterobranchia"].reset_index()
print(df_heterobranchia)

# -- -- -- --
# x
# x
# y
# w
# x
# -- -- 
# #sequencia = input("Informe a sequência a ser classificada\n")
# sequencia = "ATTAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAATGGCTACCACTTCCATGGAAGGCAGCAGGCGCGTAAATTACTCGGTGCCGATACGGCGAGGTAGTGACGACAAATACCAATGCCTACTTGCCCTTCACAGGGCGGAGGCAATTGGAATGGGCACATTGCAAATAGCTGTGCGAGTAACAATTGGAGGGAAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAG"

# k_cl = getKmersSeq(sequencia, 8)
# k_db = getKmers(organismos, 8)

# # # lendo k_mers do database PR2 (caso haja um json com o database pré-processado)
# # with open('teste.json') as f:
# #     k_db = json.load(f)


# # gerando classificao 
# classificacao = assign_family(k_cl, k_db)
# classificacao_or = assign_family_or_op(classificacao, k_cl, all_families)

# print("Classificaçao gerada")

# # escrevendo matriz familia op or
# df = (pd.DataFrame(classificacao_or)).transpose().sort_index(axis=1)
# df['results'] = df.iloc[:, 1:].sum(axis=1)
# df.sort_values(by=['results'],ascending=False, inplace=True)

# print(df['results'].head(50))

# #df.to_excel('66OCH0069500812.xlsx', index=True)


# # df2 = df.iloc[:, -1]
# # df2 = df2.head(6)
# # fig, ax = plt.subplots(figsize=(10, 5))
# # bars = ax.barh(df2.index,df2.results,color='#6892ff')
# # ax.bar_label(bars)
# # ax.invert_yaxis()
# # #ax.set_xlabel('Frequency')
# # #ax.set_ylabel('Organism')
# # ax.set_xticks([])
# # ax.set_facecolor('white')
# # plt.show()
