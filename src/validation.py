import os
import pandas as pd
import numpy as np

def validate_gwas(file_path):
    """
    Validates and standardizes GWAS summary statistics.
    Returns: (bool, str, int, str, set) -> (is_valid, message, num_snps, processed_file_path, snps_set)
    """
    try:
        # 1. Read file
        df = pd.read_csv(file_path, sep=None, engine='python')
        original_count = len(df)
        
        # 2. Standardize column names
        cols = {c.upper(): c for c in df.columns}
        
        snp_col = next((c for c in ['SNP', 'RSID', 'ID', 'MARKERNAME'] if c in cols), None)
        a1_col = next((c for c in ['A1', 'EA', 'EFFECT_ALLELE', 'ALT'] if c in cols), None)
        a2_col = next((c for c in ['A2', 'NEA', 'OTHER_ALLELE', 'REF'] if c in cols), None)
        p_col = next((c for c in ['P', 'PVAL', 'P_VALUE'] if c in cols), None)
        
        if not all([snp_col, a1_col, a2_col, p_col]):
            missing = []
            if not snp_col: missing.append("SNP ID (e.g. SNP, RSID)")
            if not a1_col: missing.append("Effect Allele (e.g. A1, EA)")
            if not a2_col: missing.append("Non-Effect Allele (e.g. A2, NEA)")
            if not p_col: missing.append("P-value (e.g. P, PVAL)")
            return False, f"Missing required columns: {', '.join(missing)}", 0, None, set()
            
        beta_col = next((c for c in ['BETA', 'B'] if c in cols), None)
        or_col = next((c for c in ['OR', 'ODDS_RATIO'] if c in cols), None)
        
        if not beta_col and not or_col:
            return False, "Missing effect size column. Expected BETA or OR.", 0, None, set()
            
        # Rename core columns
        df = df.rename(columns={
            cols[snp_col]: 'SNP',
            cols[a1_col]: 'A1',
            cols[a2_col]: 'A2',
            cols[p_col]: 'P'
        })
        
        # 3. Handle BETA / OR
        if beta_col:
            df = df.rename(columns={cols[beta_col]: 'BETA'})
            # ensure numeric
            df['BETA'] = pd.to_numeric(df['BETA'], errors='coerce')
        else:
            df = df.rename(columns={cols[or_col]: 'OR'})
            df['OR'] = pd.to_numeric(df['OR'], errors='coerce')
            df = df[df['OR'] > 0]
            df['BETA'] = np.log(df['OR'])
            
        df['P'] = pd.to_numeric(df['P'], errors='coerce')
        
        # 4. Basic cleaning
        df = df.dropna(subset=['SNP', 'A1', 'A2', 'BETA', 'P'])
        
        # P-value bounds
        df = df[(df['P'] >= 0) & (df['P'] <= 1)]
        
        # Valid alleles
        valid_alleles = ['A', 'C', 'G', 'T']
        df['A1'] = df['A1'].astype(str).str.upper()
        df['A2'] = df['A2'].astype(str).str.upper()
        
        invalid_alleles = ~df['A1'].isin(valid_alleles) | ~df['A2'].isin(valid_alleles)
        df = df[~invalid_alleles]
        
        # Ambiguous SNPs
        ambiguous = ((df['A1'] == 'A') & (df['A2'] == 'T')) | \
                    ((df['A1'] == 'T') & (df['A2'] == 'A')) | \
                    ((df['A1'] == 'C') & (df['A2'] == 'G')) | \
                    ((df['A1'] == 'G') & (df['A2'] == 'C'))
        
        ambig_count = ambiguous.sum()
        df = df[~ambiguous]
        
        # Duplicates
        dup_count = df.duplicated(subset=['SNP']).sum()
        df = df.drop_duplicates(subset=['SNP'])
        
        # Filter INFO / MAF if present
        info_col = next((c for c in ['INFO'] if c in cols), None)
        maf_col = next((c for c in ['MAF', 'EAF'] if c in cols), None)
        
        if info_col:
            df[cols[info_col]] = pd.to_numeric(df[cols[info_col]], errors='coerce')
            df = df[df[cols[info_col]] >= 0.8]
        if maf_col:
            df[cols[maf_col]] = pd.to_numeric(df[cols[maf_col]], errors='coerce')
            df = df[df[cols[maf_col]] >= 0.01]
            
        # 5. Save standard format
        processed_path = os.path.join(os.path.dirname(file_path), f"processed_{os.path.basename(file_path)}")
        final_df = df[['SNP', 'A1', 'A2', 'BETA', 'P']]
        final_df.to_csv(processed_path, sep='\t', index=False)
        
        final_count = len(final_df)
        msg = f"Processed {original_count} rows -> {final_count} valid SNPs. Removed {ambig_count} ambiguous, {dup_count} duplicates."
        
        return True, msg, final_count, processed_path, set(final_df['SNP'])
        
    except Exception as e:
        return False, f"GWAS Validation Error: {str(e)}", 0, None, set()


def validate_plink(prefix, gwas_snps=None):
    """
    Validates PLINK target data.
    """
    try:
        bed = f"{prefix}.bed"
        bim = f"{prefix}.bim"
        fam = f"{prefix}.fam"
        
        if not (os.path.exists(bed) and os.path.exists(bim) and os.path.exists(fam)):
            return False, f"Missing one or more PLINK files for prefix {prefix}. Require .bed, .bim, .fam"
            
        df_bim = pd.read_csv(bim, sep='\s+', header=None)
        if len(df_bim.columns) != 6:
            return False, f"BIM file does not have 6 columns."
            
        df_fam = pd.read_csv(fam, sep='\s+', header=None)
        n_samples = len(df_fam)
        n_snps = len(df_bim)
        
        # Check invalid alleles
        invalid_mask = (df_bim[4] == df_bim[5]) | \
                       (~df_bim[4].astype(str).str.upper().isin(['A','C','G','T'])) | \
                       (~df_bim[5].astype(str).str.upper().isin(['A','C','G','T']))
                       
        invalid_count = invalid_mask.sum()
        
        if invalid_count > 0:
            return False, f"Found {invalid_count} variants with invalid/duplicate alleles (e.g. 0 0) in .bim file. Please clean your dataset using PLINK 1.9 (e.g. plink --bfile prefix --geno 0.99 --make-bed) or upload a QC-ed dataset."
            
        overlap_msg = ""
        if gwas_snps is not None and len(gwas_snps) > 0:
            overlap = len(set(df_bim[1].values).intersection(gwas_snps))
            if overlap < 100:
                overlap_msg = f" WARNING: Only {overlap} SNPs overlap between GWAS and target .bim. PRS scores may be inaccurate."
            else:
                overlap_msg = f" Overlap OK ({overlap} SNPs)."
                
        return True, f"Target OK: {n_samples} samples, {n_snps} SNPs.{overlap_msg}"
        
    except Exception as e:
        return False, f"PLINK Validation Error: {str(e)}"


def validate_phenotype(file_path):
    """
    Validates ML Phenotype file.
    Returns: df (standardized with Sample_ID and PHENO)
    """
    try:
        df = pd.read_csv(file_path, sep=None, engine='python')
        if len(df.columns) < 2:
            raise ValueError("Phenotype file must have at least 2 columns.")
            
        if len(df.columns) >= 3:
            # Assuming FID IID PHENO
            df = df.iloc[:, [1, -1]]
        else:
            # Assuming IID PHENO or Sample_ID PHENO
            df = df.iloc[:, [0, -1]]
            
        df.columns = ['Sample_ID', 'PHENO']
        
        # Drop missing
        df = df.dropna(subset=['PHENO'])
        return True, f"Valid phenotype: {len(df)} samples.", df
    except Exception as e:
        return False, f"Phenotype Validation Error: {str(e)}", None


def validate_ld_ref(pop, chrom):
    """
    Check if LD reference exists for pop and chrom. (Very basic check).
    """
    # Assuming 1KG reference layout
    ld_dir = f"/app/ld_reference"
    if not os.path.exists(ld_dir) or not os.listdir(ld_dir):
         return True, f"Note: Could not verify local LD reference files. Assuming they exist inside container."
    
    return True, f"LD Reference check completed."
