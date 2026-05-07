import streamlit as st

def render_sidebar():
    """Renders the sidebar and returns user configurations."""
    st.sidebar.title("PRS Configuration")
    
    # 1. Mode Selector
    mode = st.sidebar.radio(
        "Operating Mode",
        options=["Single Ancestry", "Multi Ancestry"],
        help="Select Multi Ancestry to use multiple GWAS datasets from different populations."
    )
    
    st.sidebar.divider()
    
    # 2. File Uploads
    st.sidebar.subheader("Input Data")
    gwas_files = st.sidebar.file_uploader(
        "Upload GWAS Summary Statistics", 
        accept_multiple_files=True if mode == "Multi Ancestry" else False,
        help="Upload one or more GWAS summary statistics files."
    )
    
    gwas_info = []
    if gwas_files:
        gwas_list = gwas_files if isinstance(gwas_files, list) else [gwas_files]
        for idx, f in enumerate(gwas_list):
            with st.sidebar.expander(f"Details for {f.name}", expanded=True):
                pop = st.selectbox("Population", ["EUR", "AFR", "EAS", "SAS", "AMR", "Mixed"], key=f"pop_{idx}")
                n_gwas = st.number_input("Sample Size (N)", min_value=1000, value=100000, step=1000, key=f"n_{idx}")
                gwas_info.append({"file": f, "pop": pop, "n_gwas": n_gwas})

    st.sidebar.divider()
    st.sidebar.subheader("Quality Control (QC) Options")

    remove_ambiguous = st.sidebar.checkbox(
        "Remove Ambiguous SNPs (A/T, C/G)",
        value=True, # Default
        help="Remove ambiguous SNPs"
    )

    st.sidebar.divider()
    
    st.sidebar.subheader("Target Data (Testing)")
    target_data = st.sidebar.file_uploader(
        "Upload Target Genotype Data (.bed/.bim/.fam)", 
        accept_multiple_files=True,
        help="Upload PLINK binary files."
    )
    
    st.sidebar.divider()
    
    st.sidebar.subheader("Validation Data (Optional)")
    val_data = st.sidebar.file_uploader(
        "Upload Validation Genotype Data (.bed/.bim/.fam)", 
        accept_multiple_files=True,
        help="Upload PLINK binary files used to fit combining weights."
    )
    val_pheno = st.sidebar.file_uploader(
        "Upload Validation Phenotype File",
        help="Phenotype for fitting weights (and evaluating ML models)."
    )
    val_covar = st.sidebar.file_uploader(
        "Upload Validation Covariates File (Optional)",
        help="Covariates to adjust for in the regression."
    )
    
    st.sidebar.divider()
    
    # 3. Ancestry Selection
    st.sidebar.subheader("Ancestry")
    ancestry_options = ["EUR", "AFR", "EAS", "SAS", "AMR", "Mixed"]
    if mode == "Single Ancestry":
        ancestry = st.sidebar.selectbox("Select Target Ancestry", ancestry_options)
    else:
        ancestry = st.sidebar.multiselect("Select Target Ancestries", ancestry_options, default=["EUR", "AFR"])
        
    st.sidebar.divider()
    
    # 4. PRS Methods Selection
    st.sidebar.subheader("PRS Methods")
    available_methods = ["PRS-CSx", "TL-PRS", "CT-SLEB", "PROSPER", "ME-BAYES SL"]
    selected_methods = st.sidebar.multiselect(
        "Select Methods to Run",
        available_methods,
        default=["PRS-CSx", "CT-SLEB"] if mode == "Multi Ancestry" else ["PRS-CSx"]
    )
    
    # 5. Hyperparameters
    params = {}
    if selected_methods:
        with st.sidebar.expander("Advanced Hyperparameters"):
            for method in selected_methods:
                st.markdown(f"**{method}**")
                if method == "PRS-CSx":
                    params[method] = {
                        "phi": st.text_input("Global shrinkage parameter (phi)", value="1e-2", key="phi"),
                        "a": st.number_input("Parameter 'a'", value=1.0, key="a")
                    }
                elif method == "TL-PRS":
                    params[method] = {
                        "lambda": st.text_input("Lambda", value="0.1", key="lambda")
                    }
                elif method == "CT-SLEB":
                    params[method] = {
                        "r2": st.text_input("LD r2 threshold", value="0.1", key="r2"),
                        "p": st.text_input("P-value threshold", value="0.05", key="p")
                    }
                elif method == "PROSPER":
                    params[method] = {
                        "alpha": st.text_input("Penalty weight (alpha)", value="0.5", key="alpha")
                    }
                elif method == "ME-BAYES SL":
                    params[method] = {
                        "iterations": st.number_input("MCMC Iterations", value=1000, key="iters")
                    }
                    
    st.sidebar.divider()
    
    # 6. ML Predictors
    st.sidebar.subheader("ML Predictors")
    available_ml = ["SVM", "GLM", "Random Forest"]
    selected_ml = st.sidebar.multiselect("Select Models to Train", available_ml, default=["GLM"])
    
    is_binary = st.sidebar.checkbox("Is phenotype binary?", value=False)
    
    return {
        "mode": mode,
        "gwas_info": gwas_info,
        "target_data": target_data,
        "val_data": val_data,
        "val_pheno": val_pheno,
        "val_covar": val_covar,
        "ancestry": ancestry,
        "selected_methods": selected_methods,
        "params": params,
        "selected_ml": selected_ml,
        "is_binary": is_binary,
        "remove_ambiguous": remove_ambiguous,
    }
