import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Cargar tabla
df = pd.read_csv('mutaciones_consolidadas.csv')

# Heatmap: Presencia de mutaciones
pivot_presence = pd.crosstab(df['Tipo'], df['Gen'])
plt.figure(figsize=(25, 10))
sns.heatmap(pivot_presence, cmap='Reds')
plt.title('Presencia de mutaciones por tipo y gen')
plt.tight_layout()
plt.savefig('heatmap_mutaciones.png')

# Conteo de mutaciones
plt.figure(figsize=(10, 6))
df['Muestra'].value_counts().plot(kind='bar')
plt.title('Número total de mutaciones por muestra')
plt.xlabel('Muestra')
plt.ylabel('Número de mutaciones')
plt.tight_layout()
plt.savefig('mutaciones_conteo.png')

# Distribución de mutación
plt.figure(figsize=(12, 6))
pd.crosstab(df['Muestra'], df['Tipo']).plot(kind='bar', stacked=True)
plt.title('Distribución de tipos de mutación por muestra')
plt.xlabel('Muestra')
plt.ylabel('Conteo')
plt.legend(title='Tipo de mutación', bbox_to_anchor=(0.70, 0.90), loc='upper left')
plt.tight_layout()
plt.savefig('mutaciones_tipo.png')

# 4. Tabla resumen por muestra
summary_table = df.groupby(['Muestra', 'Gen']).agg({
    'Tipo': 'count',
    'AF': 'mean',
    'Impacto': lambda x: (x == 'Pathogénico').sum()
}).rename(columns={'Tipo': 'Count', 'AF': 'AF_promedio', 'Impacto': 'Pathogenic_count'})

summary_table.to_csv('resumen_por_muestra.csv')
print("Gráficos y tablas guardados.")