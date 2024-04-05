from collections import defaultdict
from functions import getKmersSeq,getKmers, assign_family, assign_family_or_op, frequency_per_family
import pandas as pd

# Open and read the database file
database = open("database/pr2_version_5.0.0_SSU_dada2.fasta", "r")
db_rows  = database.readlines()

# Initialize the datastructures
db_organisms = defaultdict(int)
all_families = set()
count = defaultdict(int)

# Parse the database rows to extract organism information
for i in range(0, len(db_rows ), 2):
  all_families.add(db_rows [i].split(";")[6])
  count[db_rows [i].split(";")[6]] +=1
  db_organisms[db_rows [i].split(";")[6] + " " + str(count[db_rows [i].split(";")[6]])] = db_rows [i+1]

# Get the sequence input from the user
sequency = input("Enter the sequence to be assigned\n")

# Extract k-mers from the input sequence and database organisms
k_cl = getKmersSeq(sequency, 8)
k_db = getKmers(db_organisms, 8)

# Assign families based on the k-mers
assignment = assign_family(k_cl, k_db)

# Assign families based on the k-mers - logical OR operation version
assignment_or = assign_family_or_op(assignment, k_cl, all_families)

print("Assignment done\n")

# Convert the assignment to a DataFrame for better visualization
df = (pd.DataFrame(assignment_or)).transpose().sort_index(axis=1)
df['Sumif 1'] = df.iloc[:, 1:].sum(axis=1)
df.sort_values(by=['Sumif 1'],ascending=False, inplace=True)

# Print the top 50 results
print(df['results'].head(50))

# Export the DataFrame to an Excel file
df.to_excel('assignment.xlsx', index=True)