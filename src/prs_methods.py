import os
import time
import pandas as pd
import numpy as np
import subprocess
import streamlit as st
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression, LogisticRegression

def _read_plink_score(file_path):
    """Reads PLINK .sscore or .profile output and returns IID and scores."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Score file not found: {file_path}")
    
    df = pd.read_csv(file_path, sep=r'\s+', engine='python')
    # plink2 typically outputs IID and SCORE1_AVG or similar
    score_col = [col for col in df.columns if 'SCORE' in col.upper() or col.upper() == 'PRS']
    if not score_col:
        # Fallback to last column
        score_col = df.columns[-1]
    else:
        score_col = score_col[0]
        
    return df[['IID', score_col]]


def _normalise_chromosomes(chromosomes):
    if chromosomes in (None, "", "1"):
        return ["1"]
    if chromosomes == "1-22":
        return [str(chrom) for chrom in range(1, 23)]
    if isinstance(chromosomes, (list, tuple)):
        return [str(chrom) for chrom in chromosomes]
    return [str(chromosomes)]


def _plink_chrom_filter(chromosomes):
    chrom_list = _normalise_chromosomes(chromosomes)
    if chrom_list == [str(chrom) for chrom in range(1, 23)]:
        return "autosome", "autosome"
    return chrom_list[0], f"chr{chrom_list[0]}"


def _prscsx_chrom_args(chromosomes):
    chrom_list = _normalise_chromosomes(chromosomes)
    if len(chrom_list) == 1:
        return [f"--chrom={chrom_list[0]}"]
    return ["--chrom"] + chrom_list

def run_prs_csx(gwas_files, gwas_pops, gwas_ns, target_plink, params, val_prefix=None, val_pheno=None, val_covar=None, is_binary=False):
    """Run PRS-CSx with robust multi-ancestry support."""
    out_dir = "/app/results/prscsx"
    
    # Clean old results to avoid stale data
    import glob
    for old_file in glob.glob(f"{out_dir}/*"):
        os.remove(old_file)
    os.makedirs(out_dir, exist_ok=True)
    out_name = "prscsx"
    
    phi = params.get('phi', '1e-2')
    a = params.get('a', '1.0')
    chromosomes = params.get("chromosomes", "1")
    chrom_list = _normalise_chromosomes(chromosomes)
    plink_chrom_filter, plink_chrom_label = _plink_chrom_filter(chromosomes)

    try:
        from src.qc import (
            clean_plink_invalid_alleles,
            find_invalid_bim_variants,
        )

        needs_allele_clean = len(find_invalid_bim_variants(target_plink)) > 0
        if not target_plink.endswith(f".prscsx_{plink_chrom_label}"):
            clean_prefix = f"{target_plink}.prscsx_{plink_chrom_label}"
            target_plink, clean_summary = clean_plink_invalid_alleles(
                target_plink,
                out_prefix=clean_prefix,
                keep_chrom=plink_chrom_filter,
            )
            st.info(
                f"Prepared target genotype for PRS-CSx chromosomes {','.join(chrom_list)}: "
                f"removed {clean_summary['num_invalid_bim_variants']} invalid BIM variants; "
                f"using `{target_plink}`."
            )
    except Exception as e:
        if hasattr(st, "warning"):
            st.warning(f"Could not pre-clean target genotype for PRS-CSx: {e}")
    
    cmd = [
        "python", "/app/tools/PRScsx/PRScsx.py",
        "--ref_dir=/app/ld_reference",
        f"--bim_prefix={target_plink}",
        f"--sst_file={','.join(gwas_files)}",
        f"--n_gwas={','.join(gwas_ns)}",
        f"--pop={','.join(gwas_pops)}",
        f"--out_dir={out_dir}",
        f"--out_name={out_name}",
        f"--phi={phi}",
        f"--a={a}",
    ] + _prscsx_chrom_args(chromosomes)
    
    # Execute with live logs
    try:
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True
        )
        
        log_placeholder = st.empty()
        full_log = ""
        
        for line in process.stdout:
            full_log += line
            # Keep only the last 20 lines to avoid UI lag
            display_log = "\n".join(full_log.splitlines()[-20:])
            log_placeholder.code(display_log)
            
        process.wait()
        if process.returncode != 0:
             raise subprocess.CalledProcessError(process.returncode, cmd, output=full_log)
             
    except Exception as e:
        st.error(f"PRS-CSx execution failed!")
        raise e
        
    # Helper to score a dataset
    def score_dataset(bfile_prefix, dataset_name):
        dfs = []
        for pop in gwas_pops:
            # Concatenate per-chromosome files into a single score file
            score_file = f"{out_dir}/{out_name}_{pop}.txt"
            chr_files = sorted([f for f in os.listdir(out_dir) if f.startswith(f"{out_name}_{pop}_pst_eff_") and f.endswith(".txt")])
            if not chr_files:
                raise FileNotFoundError(f"No PRS-CSx posterior files found for {pop}")
            
            # PRS-CSx output has NO header. Format: CHR SNP BP A1 A2 BETA (6 columns, tab-separated)
            import shutil
            with open(score_file, 'wb') as outfile:
                for chr_file in chr_files:
                    filepath = os.path.join(out_dir, chr_file)
                    with open(filepath, 'rb') as infile:
                        shutil.copyfileobj(infile, outfile)
                        
            # Fast SNP counting
            with open(score_file, 'rb') as f:
                total_snps = sum(1 for _ in f)
            
            st.info(f"Combined {len(chr_files)} chromosome files for {pop}: {total_snps} SNPs total")
            
            if total_snps == 0:
                raise ValueError(f"PRS-CSx produced 0 SNPs for {pop}. Check that GWAS SNP IDs match the LD reference panel.")
            
            # PRS-CSx headerless output columns (1-indexed for PLINK2):
            # Col 1: CHR, Col 2: SNP, Col 3: BP, Col 4: A1, Col 5: A2, Col 6: BETA
            idx_snp, idx_a1, idx_beta = 2, 4, 6

            plink_out = f"{out_dir}/{dataset_name}_{pop}_scored"
            plink_cmd = [
                "plink2",
                "--bfile", bfile_prefix,
                "--score", score_file, str(idx_snp), str(idx_a1), str(idx_beta),
                "--out", plink_out
            ]
            
            try:
                subprocess.run(plink_cmd, check=True, capture_output=True, text=True)
            except subprocess.CalledProcessError as e:
                st.error(f"PLINK2 scoring failed for {pop}! stderr:\n{e.stderr}")
                raise e
            
            sscore_path = f"{plink_out}.sscore"
            df_pop = _read_plink_score(sscore_path)
            df_pop = df_pop.rename(columns={df_pop.columns[1]: f"PRS_CSx_{pop}"})
            dfs.append(df_pop)
            
        # Merge all populations by IID
        res_df = dfs[0]
        for df_pop in dfs[1:]:
            res_df = pd.merge(res_df, df_pop, on='IID', how='outer')
        return res_df

    # Score Target
    target_scores = score_dataset(target_plink, "target")
    
    # Compute combined if validation is provided
    if val_prefix and val_pheno:
        # Validation dataset provided but it's single ancestry -> bypass regression!
        if len(gwas_pops) == 1:
            target_scores['PRS_CSx_combined'] = target_scores[f"PRS_CSx_{gwas_pops[0]}"]
            return target_scores
            
        val_scores = score_dataset(val_prefix, "val")
        
        # Read Phenotype
        pheno_df = pd.read_csv(val_pheno, sep=r'\s+|,', engine='python')
        if len(pheno_df.columns) < 2:
            raise ValueError("Validation phenotype file must have at least 2 columns: IID and Phenotype")
        # Assume IID is the first column or column named IID
        if 'IID' in pheno_df.columns:
            pheno_df = pheno_df.set_index('IID')
        else:
            pheno_df = pheno_df.set_index(pheno_df.columns[0])
            pheno_df.index.name = 'IID'
        pheno_col = pheno_df.columns[-1]
        
        # Merge validation scores with phenotype
        val_merged = pd.merge(val_scores, pheno_df[[pheno_col]], on='IID', how='inner')
        if len(val_merged) == 0:
            raise ValueError("No overlapping IIDs between validation genotype and validation phenotype!")
            
        # Covariates (Optional)
        covar_cols = []
        if val_covar:
            covar_df = pd.read_csv(val_covar, sep=r'\s+|,', engine='python')
            if 'IID' in covar_df.columns:
                covar_df = covar_df.set_index('IID')
            else:
                covar_df = covar_df.set_index(covar_df.columns[0])
                covar_df.index.name = 'IID'
            covar_cols = list(covar_df.columns)
            val_merged = pd.merge(val_merged, covar_df, on='IID', how='inner')
            
        # Standardize ancestry-specific PRS using validation set mean and SD
        prs_cols = [f"PRS_CSx_{p}" for p in gwas_pops]
        scaler = StandardScaler()
        val_merged[prs_cols] = scaler.fit_transform(val_merged[prs_cols])
        
        X = val_merged[prs_cols + covar_cols]
        y = val_merged[pheno_col]
        
        if is_binary:
            model = LogisticRegression(max_iter=1000)
        else:
            model = LinearRegression()
            
        model.fit(X, y)
        
        # Apply to target
        target_scores[prs_cols] = scaler.transform(target_scores[prs_cols])
        
        # We need to fill missing covariates in target if they were used.
        # But we only requested validation covariates. For now, assume covariates are not in target.
        # If covariates are used, predicting on target without covariates will fail unless we impute or provide target covariates.
        # Since target_covar was not in UI, we can only use the PRS weights.
        # So we extract the coefficients for the PRS columns only to form the combined score.
        prs_weights = model.coef_[0][:len(prs_cols)] if is_binary else model.coef_[:len(prs_cols)]
        target_scores['PRS_CSx_combined'] = target_scores[prs_cols].dot(prs_weights)
        
    return target_scores

def run_tl_prs(gwas_files, target_plink, params):
    raise NotImplementedError("TL-PRS mock implementation disabled in real mode.")

def run_ct_sleb(gwas_files, target_plink, params):
    raise NotImplementedError("CT-SLEB mock implementation disabled in real mode.")

def run_prosper(gwas_files, target_plink, params):
    # PROSPER has a real implementation, keeping it minimal to not break it
    out_dir = "/app/results/prosper"
    os.makedirs(out_dir, exist_ok=True)
    alpha = params.get('alpha', '0.5')
    cmd = [
        "Rscript", "/app/tools/PROSPER/PROSPER.R",
        "--gwas", gwas_files[0],
        "--target", target_plink,
        "--alpha", str(alpha),
        "--out", f"{out_dir}/prosper_out"
    ]
    subprocess.run(cmd, check=True)
    out_file = f"{out_dir}/prosper_out.txt"
    if not os.path.exists(out_file):
         raise FileNotFoundError("PROSPER failed to produce output")
    return _read_plink_score(out_file)

def run_me_bayes_sl(gwas_files, target_plink, params):
    raise NotImplementedError("ME-BAYES SL mock implementation disabled in real mode.")

def execute_prs_pipeline(methods, mode, gwas_files, gwas_pops, gwas_ns, target_file, params_dict, val_prefix=None, val_pheno=None, val_covar=None, is_binary=False):
    results_df = None
    total_methods = len(methods)
    progress_bar = st.progress(0)
    status_text = st.empty()

    for idx, method in enumerate(methods):
        status_text.text(f"Running {method} ({idx+1}/{total_methods})...")
        params = params_dict.get(method, {})
        
        try:
            if method == "PRS-CSx":
                score_df = run_prs_csx(gwas_files, gwas_pops, gwas_ns, target_file, params, val_prefix, val_pheno, val_covar, is_binary)
            elif method == "TL-PRS":
                score_df = run_tl_prs(gwas_files, target_file, params)
            elif method == "CT-SLEB":
                score_df = run_ct_sleb(gwas_files, target_file, params)
            elif method == "PROSPER":
                score_df = run_prosper(gwas_files, target_file, params)
            elif method == "ME-BAYES SL":
                score_df = run_me_bayes_sl(gwas_files, target_file, params)
            else:
                raise ValueError(f"Unknown method {method}")
                
            # Rename score column if it's PROSPER. PRS-CSx handles its own renaming
            if method != "PRS-CSx":
                score_col = [col for col in score_df.columns if col != 'IID'][0]
                score_df = score_df.rename(columns={score_col: method})
            
            if results_df is None:
                results_df = score_df
            else:
                results_df = pd.merge(results_df, score_df, on='IID', how='outer')
                
        except NotImplementedError as e:
            st.error(f"Skipping {method}: {str(e)}")
        except Exception as e:
            import traceback
            tb = traceback.format_exc()
            st.error(f"Error running {method}: {str(e)}\n\nFull traceback:\n```\n{tb}\n```")
            
        progress_bar.progress((idx + 1) / total_methods)
        
    status_text.text("All PRS methods completed successfully!")
    time.sleep(1)
    status_text.empty()
    progress_bar.empty()
    
    if results_df is not None:
        results_df = results_df.rename(columns={'IID': 'Sample_ID'})
    else:
        results_df = pd.DataFrame(columns=['Sample_ID'] + methods)
        
    return results_df
