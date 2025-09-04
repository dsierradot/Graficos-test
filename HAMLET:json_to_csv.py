import json
import pandas as pd
import glob
import os

def process_json_file(file_path):
    ##procesa y carga un json
    with open(file_path, 'r') as f:
        data = json.load(f)

    sample_name = data['metadata']['sample_name']
    mutations_list = []
    
    # Extraer mutaciones SNV/Indels
    if 'snv_indels' in data['modules'] and 'genes' in data['modules']['snv_indels']:
        snv_indels = data['modules']['snv_indels']['genes']
        for gene, variants in snv_indels.items():
            for var in variants:
                mutations_list.append({
                    'Muestra': sample_name,
                    'Gen': gene,
                    'Tipo': var.get('INFO', {}).get('TYPE', ''),
                    'Consecuencia': var.get('most_severe_consequence', ''),
                    'AF': var.get('FORMAT', {}).get('AF', ''),
                    'DP': var.get('FORMAT', {}).get('DP', ''),
                    'VD': var.get('FORMAT', {}).get('VD', ''),
                    'HGVSg': next((conseq.get('hgvsg', '') for conseq in var.get('transcript_consequences', []) if 'hgvsg' in conseq), ''),
                    'HGVSp': next((conseq.get('hgvsp', '') for conseq in var.get('transcript_consequences', []) if 'hgvsp' in conseq), ''),
                    'Impacto': 'Pathogénico' if any('clin_sig' in str(coloc) for coloc in var.get('colocated_variants', [])) else 'VUS'
                })
    
    # Extraer fusiones
    if 'fusion' in data['modules'] and 'events' in data['modules']['fusion']:
        fusions = data['modules']['fusion']['events']
        for fusion in fusions:
            if fusion.get('confidence') == 'high':  # Solo fusiones de alta confianza
                mutations_list.append({
                    'Muestra': sample_name,
                    'Gen': f"{fusion.get('gene1', '')}-{fusion.get('gene2', '')}",
                    'Tipo': 'Fusión',
                    'Consecuencia': fusion.get('reading_frame', ''),
                    'AF': 'N/A',
                    'DP': fusion.get('split_reads1', 0) + fusion.get('split_reads2', 0),
                    'VD': 'N/A',
                    'HGVSg': f"{fusion.get('breakpoint1', '')} - {fusion.get('breakpoint2', '')}",
                    'HGVSp': (fusion.get('peptide_sequence', '')[:50] + '...') if fusion.get('peptide_sequence') else 'N/A',
                    'Impacto': 'Pathogénico' if fusion.get('tags') and 'Mitelman' in str(fusion.get('tags')) else 'VUS'
                })
    
    # Extraer ITDs de FLT3
    if 'itd' in data['modules'] and 'flt3' in data['modules']['itd']:
        itd_table = data['modules']['itd']['flt3'].get('table', [])
        for itd in itd_table:
            mutations_list.append({
                'Muestra': sample_name,
                'Gen': 'FLT3',
                'Tipo': 'ITD',
                'Consecuencia': itd.get('boundary_type', ''),
                'AF': 'N/A',
                'DP': itd.get('rose_start_count', 0) + itd.get('rose_end_count', 0),
                'VD': 'N/A',
                'HGVSg': f"Start: {itd.get('td_starts', [])}, End: {itd.get('td_ends', [])}",
                'HGVSp': 'N/A',
                'Impacto': 'Pathogénico'
            })
    
    # Extraer ITDs de KMT2A si existen
    if 'itd' in data['modules'] and 'kmt2a' in data['modules']['itd']:
        kmt2a_table = data['modules']['itd']['kmt2a'].get('table', [])
        for itd in kmt2a_table:
            mutations_list.append({
                'Muestra': sample_name,
                'Gen': 'KMT2A',
                'Tipo': 'ITD',
                'Consecuencia': itd.get('boundary_type', ''),
                'AF': 'N/A',
                'DP': itd.get('rose_start_count', 0) + itd.get('rose_end_count', 0),
                'VD': 'N/A',
                'HGVSg': f"Start: {itd.get('td_starts', [])}, End: {itd.get('td_ends', [])}",
                'HGVSp': 'N/A',
                'Impacto': 'Pathogénico'
            })
    
    return mutations_list

def main():
    # Buscar todos los archivos JSON en el directorio actual
    json_files = glob.glob('*.json')
    
    if not json_files:
        print("No se encontraron archivos JSON en el directorio actual.")
        return
    
    print(f"Procesando {len(json_files)} archivos JSON...")
    
    all_mutations = []
    
    for json_file in json_files:
        try:
            print(f"Procesando: {json_file}")
            mutations = process_json_file(json_file)
            all_mutations.extend(mutations)
        except Exception as e:
            print(f"Error procesando {json_file}: {e}")
    
    # Crear DataFrame consolidado
    if all_mutations:
        df = pd.DataFrame(all_mutations)
        
        # Guardar en CSV
        output_file = 'mutaciones_consolidadas.csv'
        df.to_csv(output_file, index=False, encoding='utf-8')
        
        # Mostrar resumen
        print(f"\nResumen de mutaciones encontradas:")
        print(f"Total de mutaciones: {len(df)}")
        print(f"Muestras procesadas: {df['Muestra'].nunique()}")
        print(f"Archivos procesados: {len(json_files)}")
        print(f"\nDatos guardados en: {output_file}")
        
        # Mostrar vista previa
        print("\nVista previa de los datos:")
        print(df.head())
        
    else:
        print("No se encontraron mutaciones en los archivos procesados.")

if __name__ == "__main__":
    main()