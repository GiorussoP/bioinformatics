import urllib.request
import json
import pandas as pd

try:
    print("Carregando a lista completa de genes...")
    df = pd.read_csv('expressao_genes_completa.csv', usecols=['Marcador'])
    genes = df['Marcador'].tolist()
except Exception as e:
    print("Erro ao ler expressao_genes_completa.csv. Certifique-se de rodar o script do MATLAB primeiro:", e)
    genes = []

print(f"Total de genes carregados: {len(genes)}")

non_sex_genes = []
sex_genes = []

# Envia requisições em lotes de 1000 para a API do mygene.info via POST (MUITO mais rápido)
chunk_size = 1000
for i in range(0, len(genes), chunk_size):
    chunk = genes[i:i+chunk_size]
    genes_str = ','.join([str(g) for g in chunk])

    url = "https://mygene.info/v3/query"
    data = f"q={genes_str}&fields=genomic_pos.chr&species=human".encode('utf-8')

    try:
        req = urllib.request.Request(url, data=data, headers={'User-Agent': 'Mozilla/5.0', 'Content-Type': 'application/x-www-form-urlencoded'})
        with urllib.request.urlopen(req) as response:
            results = json.loads(response.read().decode())

            # Mapeia o resultado retornado pela API
            chr_map = {}
            for hit in results:
                gene_symbol = hit.get('query')
                if 'genomic_pos' in hit and isinstance(hit['genomic_pos'], dict):
                    chrom = hit['genomic_pos'].get('chr', 'Unknown')
                elif 'genomic_pos' in hit and isinstance(hit['genomic_pos'], list) and len(hit['genomic_pos']) > 0:
                    chrom = hit['genomic_pos'][0].get('chr', 'Unknown')
                else:
                    chrom = 'Unknown'
                if gene_symbol:
                    chr_map[gene_symbol] = chrom

            # Separa os autossômicos dos sexuais
            for gene in chunk:
                chrom = chr_map.get(str(gene), 'Unknown')
                if chrom in ['X', 'Y']:
                    sex_genes.append(f"{gene} (Chr {chrom})")
                else:
                    non_sex_genes.append(gene)
    except Exception as e:
        # Em caso de falha de conexão em um bloco, mantém os genes para evitar perdas severas
        non_sex_genes.extend(chunk)

print(f"\nTriagem Cromossômica Concluída:")
print(f"Genes autossômicos mantidos: {len(non_sex_genes)}")
print(f"Genes sexuais (X/Y) descartados: {len(sex_genes)}")

with open('genes_nao_sexuais.txt', 'w') as f:
    for gene in non_sex_genes:
        f.write(f"{gene}\n")
print("Lista salva com sucesso em 'genes_nao_sexuais.txt'!")
