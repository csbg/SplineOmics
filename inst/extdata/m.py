import pandas as pd
import random

# Load the Excel file
file_path = 'proteomics_data.xlsx'
df = pd.read_excel(file_path)

# Load the TSV file
tsv_path = 'ncbi_dataset.tsv'
tsv_df = pd.read_csv(tsv_path, sep='\t')

# Extract 'Symbol' and 'Gene name' columns from the TSV file
gene_symbols = tsv_df['Symbol'].tolist()
gene_names = tsv_df['Gene name'].tolist()

# Rename the 'Gene' column to 'Gene symbol'
df.rename(columns={'Gene': 'Gene symbol'}, inplace=True)

# Add a new 'Gene name' column initialized with NaN
df['Gene name'] = pd.NA

# Modify the 'Gene symbol' and 'Gene name' columns starting from the 4th row
for i in range(3, len(df)):
    random_index = random.randint(0, len(gene_symbols) - 1)
    df.at[i, 'Gene symbol'] = gene_symbols[random_index]
    df.at[i, 'Gene name'] = gene_names[random_index]

# Modify the 'Protein' column starting from the 4th row
for i in range(3, len(df)):
    df.at[i, 'Protein'] = f'Protein {i - 2}'

# Save the modified DataFrame back to Excel
df.to_excel(file_path, index=False)

print("File updated successfully with gene symbols and gene names.")

