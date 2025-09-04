import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Cargar tabla
df = pd.read_csv('mutaciones_consolidadas.csv')

genes_interes = {"NPM1", "KIT", "TP53", "KRAS", "WT1", "NRAS", "CEBPA", "PTPN11", "DNMT3A", "GATA2", "IDH1", "MYC", "IDH2", "MPL", "ZRSR2", "SMC3", "U2AF1", "SRSF2", "SF3B1", "ASXL1", "STAG2", "BCOR", "RUNX1", "EZH2", "MLL", "PHF6", "SF1", "NF1", "CUX1", "SETBP1", "FLT3", "TET2"}

df['Gen_clean'] = df['Gen'].str.extract(r'^([A-Za-z0-9]+)')[0]

# Filtrar solo genes de interés
df_filtrado = df[df['Gen_clean'].isin(genes_interes)]

# Heatmap: tratar de especificar los cromosomas para hacer un chr vs gen
pivot_presence = pd.crosstab(df['Muestra'], df['Gen_clean'])
plt.figure(figsize=(18, 20))
sns.heatmap(pivot_presence, cmap='Reds')
plt.title('Presencia de mutaciones por gen y tipo')
plt.tight_layout()
plt.savefig('heatmap_test.png')


# Conteo de mutaciones
plt.figure(figsize=(10, 6))
df['HGVSg'].str.split(":").str[0].value_counts().plot(kind='bar')
plt.title('Número total de mutaciones por tipo')
plt.xlabel('HGVSg')
plt.ylabel('Número de mutaciones')
plt.tight_layout()
plt.savefig('mutaciones_test.png')