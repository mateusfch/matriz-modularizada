from collections import defaultdict

# k_mers - para o database
def getKmers(organismos, tamanho):
  k_mers = {}
  for genero, sequencia in organismos.items():
    lista_k_mers = {}
    for i in range (0, len(sequencia)-tamanho+1):
      lista_k_mers[sequencia[i:i+tamanho]] = lista_k_mers.get(sequencia[i:i+tamanho],0) + 1
    k_mers[genero] = lista_k_mers
  return k_mers

# k_mers - para a sequência a ser classificada
def getKmersSeq(sequencia, tamanho):
  k_mers = {}
  for i in range(0, len(sequencia)-tamanho+1):
    k_mers[sequencia[i:i+tamanho]] = k_mers.get(sequencia[i:i+tamanho],0) + 1
  return k_mers


# classificar familia do organismo
def assign_family(k_cl, k_db):
  aparece = defaultdict(defaultdict)
  for k_dbKEY, k_dbVAL in k_db.items():
    intersec = k_dbVAL.keys() & k_cl.keys()
    aparece_KMER = defaultdict(int)
    for k_clKEY in k_cl.keys():
      aparece_KMER[k_clKEY] = 1 if k_clKEY in intersec else 0 
    aparece[k_dbKEY] = aparece_KMER
  return aparece

# classificar familia do organismo usando versao or
def assign_family_or_op(aparece, k_cl, all_families):
  aparece2 = {}
  for fam_nome in all_families:
    aparece2[fam_nome] = {k_mer: 0 for k_mer in k_cl.keys()}
  for a_key, a_val in aparece.items():
    fam_nome = a_key.split(" ")[0]
    for k_mer, freq in a_val.items():
      aparece2[fam_nome][k_mer] = 1 if freq == 1 else 0
  return aparece2