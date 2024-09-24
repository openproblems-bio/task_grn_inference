import anndata as ad 
import pandas as pd
import numpy as np
import scanpy as sc
import biomart


## VIASH START
par = {
    'gene_annotations': 'resources/prior/gene_annotations.txt'
}
## VIASH END

server = biomart.BiomartServer("http://www.ensembl.org/biomart")

# Use the Ensembl Genes mart
ensembl = server.datasets['hsapiens_gene_ensembl']

# List available attributes and filters
print(ensembl.attributes)
print(ensembl.filters)

# Define the attributes 
attributes = [
    'ensembl_gene_id',
    'external_gene_name',
    'gene_biotype',
    'chromosome_name',
    'start_position',
    'end_position',
    'strand'
]

# Query BioMart
results = ensembl.search({
    'filters': {},
    'attributes': attributes
})

results_content = results.content.decode('utf-8').splitlines()

save the results to a file
with open(par['gene_annotations'], 'w') as f:
    for line in results_content:
        f.write(line + '\n')
