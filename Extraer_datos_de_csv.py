import pandas as pd

# Cargar tabla #
df = pd.read_csv('mutaciones_consolidadas.csv')
muestra = "L0243"

genes_interes = {"NPM1", "KIT", "TP53", "KRAS", "WT1", "NRAS", "CEBPA", "PTPN11", "DNMT3A", "GATA2", "IDH1", "MYC", "IDH2", "MPL", "ZRSR2", "SMC3", "U2AF1", "SRSF2", "SF3B1", "ASXL1", "STAG2", "BCOR", "RUNX1", "EZH2", "MLL", "PHF6", "SF1", "NF1", "CUX1", "SETBP1", "FLT3", "TET2"}

# busca genes segun 3 posibles mutaciones presentes #
genes_mutados = df[(df['Muestra'] == muestra) & (df['Gen'].isin(genes_interes))]['Gen'].unique().tolist()
genes_fusionados = df[(df['Muestra'] == muestra) & (df['Tipo'].str.contains('Fusión'))]['HGVSg'].unique().tolist()
itd = df[(df['Muestra'] == muestra) & (df['Tipo'].str.contains('ITD'))]['Gen'].unique().tolist()

print("")
print(f"Muestra: {muestra}")
print(f"Mutados: {genes_mutados}")
print(f"Fusiones: {genes_fusionados}")
print(f"ITD: {itd}")