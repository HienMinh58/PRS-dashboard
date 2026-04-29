import pandas as pd
import numpy as np

print("Loading reference panel...")
ref_df = pd.read_csv(r"d:\prs\ld_reference\snpinfo_mult_1kg_hm3", sep=r'\s+', nrows=100000)

print("Fixing TOY_BASE_GWAS.assoc...")
gwas = pd.read_csv(r"d:\prs\data\TOY_BASE_GWAS.assoc", sep=r'\s+')
n_gwas = len(gwas)

# Map real rsIDs and alleles from reference
gwas['SNP'] = ref_df['SNP'].iloc[:n_gwas].values
gwas['A1'] = ref_df['A1'].iloc[:n_gwas].values
gwas['A2'] = ref_df['A2'].iloc[:n_gwas].values

# Convert OR to BETA (log-odds) for consistency
gwas['BETA'] = np.log(gwas['OR'].replace(0, 0.001))  # avoid log(0)

# PRS-CSx expects columns in EXACT order: SNP A1 A2 BETA P
# No other columns should be present!
gwas_out = gwas[['SNP', 'A1', 'A2', 'BETA', 'P']]
gwas_out.to_csv(r"d:\prs\data\TOY_BASE_GWAS.assoc", sep='\t', index=False)
print(f"GWAS fixed: {len(gwas_out)} SNPs, columns = {list(gwas_out.columns)}")
print(gwas_out.head())

print("\nFixing TOY_TARGET_DATA.bim...")
bim = pd.read_csv(r"d:\prs\data\TOY_TARGET_DATA.bim", sep=r'\s+', header=None)
n_bim = len(bim)
bim[1] = ref_df['SNP'].iloc[:n_bim].values
bim[0] = ref_df['CHR'].iloc[:n_bim].values
bim[3] = ref_df['BP'].iloc[:n_bim].values
bim[4] = ref_df['A1'].iloc[:n_bim].values
bim[5] = ref_df['A2'].iloc[:n_bim].values
bim.to_csv(r"d:\prs\data\TOY_TARGET_DATA.bim", sep='\t', header=False, index=False)

print("Fixing val_TOY_TARGET_DATA.bim...")
bim.to_csv(r"d:\prs\data\val_TOY_TARGET_DATA.bim", sep='\t', header=False, index=False)

print("\nDone! Verifying output...")
check = pd.read_csv(r"d:\prs\data\TOY_BASE_GWAS.assoc", sep=r'\s+', nrows=5)
print(check)
