from collections import defaultdict
import pandas as pd

def getKmers(organismos, tamanho):
  """
  Extract k-mers from sequences of organisms.

  Parameters:
  organisms (dict): Dictionary where keys are organism names and values are sequences.
  size (int): Size of k-mers to extract.

  Returns:
  dict: A dictionary where keys are organism names and values are dictionaries of k-mer counts.
  """
  
  k_mers = {}
  for genero, sequencia in organismos.items():
    lista_k_mers = {}
    for i in range (0, len(sequencia)-tamanho+1):
      lista_k_mers[sequencia[i:i+tamanho]] = lista_k_mers.get(sequencia[i:i+tamanho],0) + 1
    k_mers[genero] = lista_k_mers
  return k_mers

def getKmersSeq(sequencia, tamanho):
  """
  Extract k-mers from a single sequence.

  Parameters:
  sequence (str): The sequence to extract k-mers from.
  size (int): Size of k-mers to extract.

  Returns:
  dict: A dictionary of k-mer counts.
  """  

  k_mers = {}
  for i in range(0, len(sequencia)-tamanho+1):
    k_mers[sequencia[i:i+tamanho]] = k_mers.get(sequencia[i:i+tamanho],0) + 1
  return k_mers

def assign_family(k_cl, k_db):
  """
  Assign families to organisms based on k-mer intersections.

  Parameters:
  k_cl (dict): Dictionary of k-mers for the organism to classify.
  k_db (dict): Dictionary of k-mers for the database organisms.

  Returns:
  dict: A dictionary where keys are database organism names and values are dictionaries
  indicating presence (1) or absence (0) of shared k-mers with the organism to classify.
  """ 
  
  aparece = defaultdict(defaultdict)
  for k_dbKEY, k_dbVAL in k_db.items():
    intersec = k_dbVAL.keys() & k_cl.keys()
    aparece_KMER = defaultdict(int)
    for k_clKEY in k_cl.keys():
      aparece_KMER[k_clKEY] = 1 if k_clKEY in intersec else 0 
    aparece[k_dbKEY] = aparece_KMER
  return aparece

def assign_family_or_op(aparece, k_cl, all_families):
  """
  Assign families to organisms based on shared k-mers using logical OR operation.

  Parameters:
  appearances (dict): Dictionary of k-mer presence/absence for database organisms.
  k_cl (dict): Dictionary of k-mers for the organism to classify.
  all_families (list): List of all possible families.

  Returns:
  dict: A dictionary where keys are family names and values are dictionaries
  indicating presence (1) or absence (0) of shared k-mers with the organism to classify.
  """
  aparece2 = {}
  for fam_nome in all_families:
    aparece2[fam_nome] = {k_mer: 0 for k_mer in k_cl.keys()}
  for a_key, a_val in aparece.items():
    fam_nome = a_key.split()[0].strip()
    for k_mer, freq in a_val.items():
      if freq == 1:
        aparece2[fam_nome][k_mer] = 1
  return aparece2


def frequency_per_family(db_rows, family):
  """
  Filter database rows by family and return as a DataFrame.

  Parameters:
  db_rows (list): List of database rows, alternating between metadata and sequences.
  family (str): The family to filter by.

  Returns:
  pandas.DataFrame: DataFrame containing rows corresponding to the specified family.
  """  
  
  dict = defaultdict(list)
  for i in range(0, len(db_rows), 2):
    dict["Domain"].append(db_rows[i].split(";")[0])
    dict["Supergroup"].append(db_rows[i].split(";")[1])
    dict["Division"].append(db_rows[i].split(";")[2])
    dict["Subdivision"].append(db_rows[i].split(";")[3])
    dict["Class"].append(db_rows[i].split(";")[4])
    dict["Order"].append(db_rows[i].split(";")[5])
    dict["Family"].append(db_rows[i].split(";")[6])
    dict["Genus"].append(db_rows[i].split(";")[7])
    dict["Species"].append(db_rows[i].split(";")[8])
    dict["Sequence"].append(db_rows[i+1])
  df = pd.DataFrame(dict)
  return df.loc[df['Family']==family].reset_index(drop=True)