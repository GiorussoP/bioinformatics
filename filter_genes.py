import urllib.request
import json
import pandas as pd

# Load genes from CSV. We are going to remove sexual dimorphysm genes since we know the variance in data could be explainmed by sex differences in expression
try:
    df = pd.read_csv('expressao_genes_reduzidos.csv')
    genes = df.iloc[:, 0].tolist()
except Exception as e:
    print("Error reading CSV:", e)
    genes = []

print(f"Total genes loaded: {len(genes)}")

non_sex_genes = []
sex_genes = []

for gene in genes:
    try:
        url = f"https://mygene.info/v3/query?q=symbol:{gene}&fields=symbol,genomic_pos.chr&species=human"
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read().decode())
            if data.get('hits'):
                chrom = data['hits'][0].get('genomic_pos', {}).get('chr', 'Unknown') if isinstance(data['hits'][0].get('genomic_pos'), dict) else 'Unknown/Multiple'
                if chrom in ['X', 'Y']:
                    sex_genes.append(f"{gene} (Chr {chrom})")
                else:
                    non_sex_genes.append(gene)
            else:
                non_sex_genes.append(gene) # Assume non-sex if not found
    except Exception as e:
        non_sex_genes.append(gene)

# Save to file so other scripts can use them dynamically
with open('genes_nao_sexuais.txt', 'w') as f:
    for g in non_sex_genes:
        f.write(g + '\n')

print(f"Foram salvos {len(non_sex_genes)} genes não-sexuais em 'genes_nao_sexuais.txt'.")
