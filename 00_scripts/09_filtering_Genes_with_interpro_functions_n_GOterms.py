import pandas as pd

# Read the TSV file
interproscan_df = pd.read_csv("/mmfs1/scratch/jatinder.singh/RNAseq_aumb_parents/13_interproscan_on_braker3_isoseqRNAseq_again/00_finally/02_all_appl/aumb_brakerRNAref_IsoSeqcomp_functional_renamed_astrix_removed_all.tsv", 
                              sep="\t", header = None, names=["protein_id", "seqMD5", "seqLength", "analysis", "signature_accession", "signature_description", "start", "stop", "score", "status", "date", "interpro_accession", "interpro_description", "go", "pathways"])

# Create a new column 'Gene' by removing the last '.<number>' part from 'protein_id'
interproscan_df["Gene"] = interproscan_df["protein_id"].str.replace(r'\.\d+$', '', regex=True)

columns_to_get = ["Gene", "interpro_description", "go"]

# Convert NaN values to empty strings in the relevant columns
introproscan_modified_df = interproscan_df[columns_to_get].copy()  # Create a copy to avoid warnings
introproscan_modified_df.loc[:, ['interpro_description', 'go']] = introproscan_modified_df[['interpro_description', 'go']].fillna('')

# Group by "Gene" and aggregate
filtered_df = introproscan_modified_df.groupby('Gene').agg({
    'interpro_description': lambda x: ', '.join(sorted(set(desc for desc in x if desc != '-'))),
    'go': lambda x: '|'.join(sorted(set(ann for ann in x if ann not in ['', '-'])))
}).reset_index()

# Filter the empty cells in interpro_description and go
filtered_df = filtered_df[(filtered_df['interpro_description'] != "") & (filtered_df['go'] != "")]
filtered_df['go'] = filtered_df['go'].str.replace('|', ', ')

# Remove duplicates within each cell of the 'go' column
filtered_df['go'] = filtered_df['go'].apply(lambda x: ', '.join(sorted(set(x.split(', ')))))
filtered_df['interpro_description'] = filtered_df['interpro_description'].apply(lambda x: ', '.join(sorted(set(x.split(', ')))))

filtered_df.to_csv("../aumb_filtered_interproscan_with_gene_function_description.csv", index=None)
