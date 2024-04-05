from collections import defaultdict
from functions import getKmersSeq,getKmers, assign_family, assign_family_or_op, frequency_per_family

database = open("database/pr2_version_5.0.0_SSU_dada2.fasta", "r")
db_rows  = database.readlines()

organismos = defaultdict(int)
count = defaultdict(int)
all_families = set()

dict = {}

for i in range(0, len(db_rows ), 2):
  all_families.add(db_rows [i].split(";")[6])
  count[db_rows [i].split(";")[6]] +=1
  organismos[db_rows [i].split(";")[6] + " " + str(count[db_rows [i].split(";")[6]])] = db_rows [i+1]

sequencia = input("Informe a sequência a ser classificada\n")

k_cl = getKmersSeq(sequencia, 8)
k_db = getKmers(organismos, 8)

classificacao = assign_family(k_cl, k_db)
classificacao_or = assign_family_or_op(classificacao, k_cl, all_families)

print("Assignment done!")

# # writing OrOp Family matrix to Excel
# df = (pd.DataFrame(classificacao_or)).transpose().sort_index(axis=1)
# df['results'] = df.iloc[:, 1:].sum(axis=1)
# df.sort_values(by=['results'],ascending=False, inplace=True)

# print(df['results'].head(50))

#df.to_excel('assignment.xlsx', index=True)