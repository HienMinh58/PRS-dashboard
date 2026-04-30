import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.figure_factory as ff

# Internal imports
from src.ui import render_sidebar
from src.prs_methods import execute_prs_pipeline
from src.ml_models import train_ml_models
from src.utils import generate_mock_phenotype, export_to_csv, export_to_excel, calculate_metrics

# --- Page Configuration ---
st.set_page_config(
    page_title="PRS Dashboard",
    page_icon=None,
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- Custom CSS: Clean light theme ---
st.markdown("""
<style>
    /* --- Import professional font --- */
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');

    /* --- Global App Backgrounds & Text --- */
    html, body, [class*="css"], .stApp {
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
        background-color: #ffffff !important;
        color: #333333 !important;
    }
    
    [data-testid="stAppViewContainer"], .main .block-container {
        background-color: #ffffff !important;
        color: #333333 !important;
    }
    
    /* --- Top Header & Toolbar --- */
    [data-testid="stHeader"] {
        background-color: #ffffff !important;
    }
    
    /* Make toolbar icons and text dark */
    [data-testid="stToolbar"] button,
    [data-testid="stStatusWidget"] label {
        color: #333333 !important;
    }
    [data-testid="stToolbar"] svg,
    [data-testid="stStatusWidget"] svg {
        fill: #333333 !important;
        stroke: #333333 !important;
    }
    
    /* Remove the colorful top decoration line */
    [data-testid="stDecoration"] {
        background-image: none !important;
        background-color: transparent !important;
    }
    
    .main .block-container {
        padding-top: 1.5rem;
        padding-bottom: 2rem;
        max-width: 1200px;
    }

    /* --- Sidebar --- */
    [data-testid="stSidebar"], [data-testid="stSidebarContent"] {
        background-color: #f7f7f7 !important;
        border-right: 1px solid #e0e0e0 !important;
    }
    [data-testid="stSidebar"] [data-testid="stMarkdownContainer"] p,
    [data-testid="stSidebar"] label,
    [data-testid="stSidebar"] span {
        font-size: 0.88rem;
        color: #444444 !important;
    }

    /* --- Headers & Markdown --- */
    h1, h2, h3, h4, h5, h6, [data-testid="stMarkdownContainer"] {
        color: #1a1a1a !important;
    }
    h1 {
        font-size: 1.6rem !important;
        font-weight: 700 !important;
        padding-bottom: 0.3rem !important;
        border-bottom: 2px solid #2c7fb8 !important;
        margin-bottom: 1rem !important;
    }
    h2 {
        font-size: 1.25rem !important;
        font-weight: 600 !important;
        margin-top: 0.8rem !important;
        margin-bottom: 0.5rem !important;
    }
    h3 {
        font-size: 1.05rem !important;
        font-weight: 600 !important;
    }
    
    p {
        color: #333333 !important;
    }

    /* --- Primary buttons: muted blue --- */
    .stButton > button[kind="primary"],
    .stButton > button[data-testid="stBaseButton-primary"] {
        background-color: #2c7fb8 !important;
        color: #ffffff !important;
        border: none !important;
        border-radius: 4px !important;
        font-weight: 500 !important;
        font-size: 0.88rem !important;
        padding: 0.45rem 1.2rem !important;
        transition: background-color 0.2s ease !important;
    }
    .stButton > button[kind="primary"]:hover,
    .stButton > button[data-testid="stBaseButton-primary"]:hover {
        background-color: #22628f !important;
    }

    /* --- Secondary buttons --- */
    .stButton > button:not([kind="primary"]):not([data-testid="stBaseButton-primary"]),
    .stDownloadButton > button {
        background-color: #ffffff !important;
        color: #333333 !important;
        border: 1px solid #cccccc !important;
        border-radius: 4px !important;
        font-size: 0.85rem !important;
    }
    .stButton > button:not([kind="primary"]):hover,
    .stDownloadButton > button:hover {
        background-color: #f0f0f0 !important;
        border-color: #999999 !important;
    }

    /* --- File uploader --- */
    [data-testid="stFileUploader"] {
        background-color: #ffffff !important;
        border: 1px dashed #b0b0b0 !important;
        border-radius: 4px !important;
        padding: 0.5rem !important;
    }
    [data-testid="stFileUploader"] section {
        background-color: #fcfcfc !important;
    }
    [data-testid="stFileUploader"] button {
        background-color: #f0f0f0 !important;
        color: #2c7fb8 !important;
        border: 1px solid #cccccc !important;
    }
    [data-testid="stFileUploader"] button:hover {
        background-color: #e0e0e0 !important;
    }
    [data-testid="stFileUploader"] small {
        color: #666666 !important;
    }

    /* --- Inputs: Text, Selectbox, Multiselect, Checkbox, Radio --- */
    input, select, textarea, 
    [data-testid="stSelectbox"] div[data-baseweb="select"],
    [data-testid="stMultiSelect"] div[data-baseweb="select"] {
        background-color: #ffffff !important;
        color: #333333 !important;
        border-color: #cccccc !important;
        border-radius: 4px !important;
    }
    
    /* Multiselect tags */
    [data-baseweb="tag"] {
        background-color: #e6f0f9 !important;
        color: #2c7fb8 !important;
        border: 1px solid #cce0f0 !important;
    }
    [data-baseweb="tag"] span {
        color: #2c7fb8 !important;
    }

    /* Radio buttons & Checkboxes */
    [data-baseweb="radio"] div[data-testid="stMarkdownContainer"] p,
    [data-baseweb="checkbox"] div[data-testid="stMarkdownContainer"] p {
        color: #333333 !important;
    }
    /* Checked state accents (overriding red/pink) */
    div[data-baseweb="radio"] > div:first-child > div[data-checked="true"] > div,
    div[data-baseweb="checkbox"] > div > div[data-checked="true"] {
        background-color: #2c7fb8 !important;
        border-color: #2c7fb8 !important;
    }

    /* --- Tabs --- */
    .stTabs [data-baseweb="tab-list"] {
        gap: 0px;
        border-bottom: 1px solid #ddd;
        background-color: transparent !important;
    }
    .stTabs [data-baseweb="tab"] {
        background-color: transparent !important;
        border: none !important;
        border-bottom: 2px solid transparent !important;
        color: #666666 !important;
        font-size: 0.88rem !important;
        font-weight: 500 !important;
        padding: 0.5rem 1rem !important;
    }
    .stTabs [aria-selected="true"] {
        border-bottom: 2px solid #2c7fb8 !important;
        color: #1a1a1a !important;
    }

    /* --- Alerts / Info boxes --- */
    [data-testid="stAlert"] {
        border-radius: 4px !important;
        font-size: 0.88rem !important;
        border-left-width: 4px !important;
        background-color: #f8f9fa !important;
        color: #333333 !important;
    }

    /* --- Dataframe --- */
    [data-testid="stDataFrame"] {
        border: 1px solid #e0e0e0 !important;
        border-radius: 4px !important;
        background-color: #ffffff !important;
    }

    /* --- Divider --- */
    hr {
        border-top: 1px solid #e8e8e8 !important;
        margin: 1rem 0 !important;
    }

    /* --- Expander --- */
    [data-testid="stExpander"] {
        border: 1px solid #e0e0e0 !important;
        border-radius: 4px !important;
        background-color: #ffffff !important;
    }
    [data-testid="stExpander"] summary p {
        color: #333333 !important;
        font-weight: 500 !important;
    }

    /* --- Caption text --- */
    .stCaption, [data-testid="stCaptionContainer"] {
        color: #777777 !important;
        font-size: 0.8rem !important;
    }
</style>
""", unsafe_allow_html=True)


# --- Main Application Logic ---
def main():
    st.title("PRS Dashboard — Multi & Single Ancestry + ML Predictors")
    st.caption("A platform for polygenic risk score calculation, multi-ancestry integration, and machine learning prediction.")
    
    # Render Sidebar and get configs
    config = render_sidebar()
    
    # Initialize Session State
    if 'prs_results' not in st.session_state:
        st.session_state.prs_results = None
    if 'ml_results' not in st.session_state:
        st.session_state.ml_results = None
    if 'phenotype' not in st.session_state:
        st.session_state.phenotype = None
    if 'ml_prs_df' not in st.session_state:
        st.session_state.ml_prs_df = None
    if 'ml_predictions' not in st.session_state:
        st.session_state.ml_predictions = None

    # ================================================================
    # SECTION 1: PRS Calculation
    # ================================================================
    st.header("PRS Calculation")
    
    col_prs_btn, col_prs_status = st.columns([1, 3])
    with col_prs_btn:
        run_prs_btn = st.button("Run PRS Calculation", type="primary", use_container_width=True)
    with col_prs_status:
        if st.session_state.prs_results is not None:
            n_samples = len(st.session_state.prs_results)
            n_scores = len([c for c in st.session_state.prs_results.columns if c != 'Sample_ID'])
            st.success(f"PRS results available: {n_samples} samples, {n_scores} score column(s)")
        else:
            st.info("No PRS results yet. Upload data and click Run PRS Calculation.")

    if run_prs_btn:
        if not config["selected_methods"]:
            st.warning("Please select at least one PRS method from the sidebar.")
        else:
            with st.spinner("Processing Genotypes and Calculating PRS..."):
                import os
                
                data_dir = "/app/data" if os.path.exists("/app") else "./data"
                os.makedirs(data_dir, exist_ok=True)
                
                from src.validation import validate_gwas, validate_plink, validate_ld_ref
                
                # --- Validate and save GWAS files ---
                gwas_paths = []
                gwas_pops = []
                gwas_ns = []
                gwas_snps_all = set()
                
                st.subheader("Data Validation")
                validation_passed = True
                
                if config["gwas_info"] and len(config["gwas_info"]) > 0:
                    for info in config["gwas_info"]:
                        f = info["file"]
                        path = os.path.join(data_dir, f.name)
                        with open(path, "wb") as out_f:
                            out_f.write(f.getbuffer())
                        
                        with st.spinner(f"Validating GWAS: {f.name}..."):
                            is_valid, msg, count, processed_path, snps = validate_gwas(path)
                        
                        if is_valid:
                            st.success(f"GWAS '{f.name}' OK: {msg}")
                            gwas_paths.append(processed_path)
                            gwas_pops.append(info["pop"])
                            gwas_ns.append(str(info["n_gwas"]))
                            gwas_snps_all.update(snps)
                        else:
                            st.error(f"GWAS '{f.name}' Failed: {msg}")
                            validation_passed = False
                else:
                    st.error("Please upload at least one GWAS summary statistics file.")
                    st.stop()
                
                # --- Validate and save target genotype ---
                target_prefix = "mock_target"
                if config["target_data"] is not None and config["target_data"] != []:
                    exts = [f.name.split('.')[-1] for f in config["target_data"]]
                    if not all(ext in exts for ext in ['bed', 'bim', 'fam']):
                        st.error("Target Data: Please upload all three PLINK binary files: .bed, .bim, .fam")
                        validation_passed = False
                    else:
                        for f in config["target_data"]:
                            path = os.path.join(data_dir, f.name)
                            with open(path, "wb") as out_f:
                                out_f.write(f.getbuffer())
                        
                        bed_file = [f.name for f in config["target_data"] if f.name.endswith('.bed')][0]
                        target_prefix = os.path.join(data_dir, bed_file.rsplit('.bed', 1)[0])
                        
                        with st.spinner("Validating Target Data..."):
                            is_valid, msg = validate_plink(target_prefix, gwas_snps_all)
                        if is_valid:
                            st.success(msg)
                        else:
                            st.error(msg)
                            validation_passed = False
                else:
                    st.error("Please upload Target Genotype Data (.bed/.bim/.fam).")
                    validation_passed = False
                
                # --- Validation data (optional) ---
                val_prefix = None
                if config["val_data"] is not None and config["val_data"] != []:
                    exts = [f.name.split('.')[-1] for f in config["val_data"]]
                    if not all(ext in exts for ext in ['bed', 'bim', 'fam']):
                        st.error("Validation Data: Please upload all three PLINK binary files for validation: .bed, .bim, .fam")
                        validation_passed = False
                    else:
                        for f in config["val_data"]:
                            path = os.path.join(data_dir, "val_" + f.name)
                            with open(path, "wb") as out_f:
                                out_f.write(f.getbuffer())
                        bed_file = [f.name for f in config["val_data"] if f.name.endswith('.bed')][0]
                        val_prefix = os.path.join(data_dir, "val_" + bed_file.rsplit('.bed', 1)[0])
                        
                        with st.spinner("Validating Validation Data..."):
                            is_valid, msg = validate_plink(val_prefix, gwas_snps_all)
                        if is_valid:
                            st.success(f"Validation Target OK: {msg}")
                        else:
                            st.error(f"Validation Target Failed: {msg}")
                            validation_passed = False
                    
                val_pheno_path = None
                if config["val_pheno"]:
                    val_pheno_path = os.path.join(data_dir, config["val_pheno"].name)
                    with open(val_pheno_path, "wb") as out_f:
                        out_f.write(config["val_pheno"].getbuffer())
                        
                val_covar_path = None
                if config["val_covar"]:
                    val_covar_path = os.path.join(data_dir, config["val_covar"].name)
                    with open(val_covar_path, "wb") as out_f:
                        out_f.write(config["val_covar"].getbuffer())

                # Check LD Reference
                if validation_passed:
                     for pop in gwas_pops:
                         valid, msg = validate_ld_ref(pop, "1")
                         if not valid:
                             st.warning(msg)
                         else:
                             st.success(f"LD Ref {pop}: {msg}")

                if not validation_passed:
                    st.error("Input validation failed. Please fix the errors above and try again.")
                    st.stop()

                # --- Execute PRS pipeline ---
                prs_df = execute_prs_pipeline(
                    methods=config["selected_methods"],
                    mode=config["mode"],
                    gwas_files=gwas_paths,
                    gwas_pops=gwas_pops,
                    gwas_ns=gwas_ns,
                    target_file=target_prefix,
                    params_dict=config["params"],
                    val_prefix=val_prefix,
                    val_pheno=val_pheno_path,
                    val_covar=val_covar_path,
                    is_binary=config["is_binary"]
                )
                
                if prs_df is None or len(prs_df) == 0:
                    st.error("No PRS results were generated. Please check the tool errors above.")
                    st.stop()
                
                if 'Sample_ID' not in prs_df.columns:
                    prs_df.insert(0, 'Sample_ID', [f"IID_{i}" for i in range(len(prs_df))])
                
                # Store PRS results
                st.session_state.prs_results = prs_df
                # Reset ML results since PRS data changed
                st.session_state.ml_results = None
                st.session_state.phenotype = None
                st.session_state.ml_prs_df = None
                
                # Auto-save PRS results to disk for later re-use
                results_dir = "/app/results" if os.path.exists("/app") else "./results"
                os.makedirs(results_dir, exist_ok=True)
                prs_save_path = os.path.join(results_dir, "prs_results.csv")
                prs_df.to_csv(prs_save_path, index=False)
                
            st.success(f"PRS Calculation Complete! {len(prs_df)} samples scored. Results saved to `{prs_save_path}`")
            st.rerun()

    # ================================================================
    # SECTION 2: ML Prediction
    # ================================================================
    st.divider()
    st.header("ML Prediction")
    
    # --- PRS Data Source ---
    st.subheader("PRS Data Source")
    prs_source = st.radio(
        "Select PRS data source:",
        options=["From current session", "Upload PRS results CSV"],
        horizontal=True,
        help="Use PRS scores from the calculation above, or upload a previously saved PRS results file."
    )
    
    ml_prs_df = None
    
    if prs_source == "From current session":
        if st.session_state.prs_results is not None:
            ml_prs_df = st.session_state.prs_results.copy()
            n_samples = len(ml_prs_df)
            n_cols = len([c for c in ml_prs_df.columns if c != 'Sample_ID'])
            st.success(f"Using session PRS results: {n_samples} samples, {n_cols} score column(s)")
        else:
            st.warning("No PRS results in current session. Please run PRS Calculation first or upload a PRS results CSV.")
    else:
        prs_csv_file = st.file_uploader(
            "Upload PRS Results CSV",
            type=["csv", "tsv", "txt"],
            help="CSV file with columns: Sample_ID (or IID) and one or more PRS score columns (e.g., PRS_CSx_EUR, PRS_CSx_AFR, PRS_CSx_combined).",
            key="ml_prs_csv_upload"
        )
        if prs_csv_file is not None:
            try:
                ml_prs_df = pd.read_csv(prs_csv_file, sep=None, engine='python')
                
                # Normalize ID column name
                if 'IID' in ml_prs_df.columns and 'Sample_ID' not in ml_prs_df.columns:
                    ml_prs_df = ml_prs_df.rename(columns={'IID': 'Sample_ID'})
                if '#FID' in ml_prs_df.columns:
                    ml_prs_df = ml_prs_df.drop(columns=['#FID'], errors='ignore')
                
                # Drop non-score columns that are common in PLINK output
                drop_cols = ['FID', 'PHENO1', 'ALLELE_CT', 'NAMED_ALLELE_DOSAGE_SUM', 'SCORE1_AVG']
                for col in drop_cols:
                    if col in ml_prs_df.columns and col != 'Sample_ID':
                        ml_prs_df = ml_prs_df.drop(columns=[col], errors='ignore')
                
                if 'Sample_ID' not in ml_prs_df.columns:
                    st.error("PRS file must contain a 'Sample_ID' or 'IID' column.")
                    ml_prs_df = None
                else:
                    st.success(f"Loaded PRS CSV: {len(ml_prs_df)} samples, columns: {list(ml_prs_df.columns)}")
            except Exception as e:
                st.error(f"Error reading PRS CSV: {e}")
                ml_prs_df = None
    
    if ml_prs_df is not None:
        st.divider()
        
        # --- Phenotype Upload ---
        st.subheader("Phenotype Data")
        ml_col1, ml_col2 = st.columns(2)
        with ml_col1:
            ml_pheno_file = st.file_uploader(
                "Upload Phenotype File",
                type=["csv", "tsv", "txt", "pheno"],
                help="File with columns: FID IID Phenotype (or IID Phenotype). Required for ML evaluation.",
                key="ml_pheno_upload"
            )
        with ml_col2:
            ml_is_binary = st.checkbox(
                "Is phenotype binary (case/control)?",
                value=config["is_binary"],
                key="ml_binary_check"
            )
            use_demo_mode = st.checkbox(
                "Demo mode (use random phenotype)",
                value=False,
                help="Generates a random phenotype for testing only. Results will NOT be scientifically meaningful."
            )
        
        # --- Feature Column Selection ---
        st.subheader("Feature Selection")
        all_score_cols = [c for c in ml_prs_df.columns if c not in ['Sample_ID', 'IID', 'FID']]
        
        if len(all_score_cols) == 0:
            st.error("No PRS score columns found in the data. Ensure your file has numeric score columns.")
        else:
            selected_features = st.multiselect(
                "Select PRS columns to use as ML features:",
                options=all_score_cols,
                default=all_score_cols,
                help="Choose which PRS score columns to include as input features for ML training."
            )
            
            if not selected_features:
                st.warning("Please select at least one PRS feature column.")
            
            st.caption(f"Selected {len(selected_features)} of {len(all_score_cols)} available feature(s)")
            
            # --- ML Model Selection ---
            selected_ml = config["selected_ml"]
            
            # --- Run ML Button ---
            run_ml_btn = st.button("Run ML Prediction", type="primary", use_container_width=False)
            
            if run_ml_btn:
                if not selected_features:
                    st.error("Please select at least one PRS feature column.")
                    st.stop()
                if not selected_ml:
                    st.error("No ML models selected. Please select at least one model from the sidebar.")
                    st.stop()
                
                prs_df = ml_prs_df.copy()
                phenotype = None
                
                # --- Load and merge phenotype ---
                if ml_pheno_file is not None:
                    from src.validation import validate_phenotype
                    ml_pheno_file.seek(0)
                    is_valid, msg, pheno_df = validate_phenotype(ml_pheno_file)
                    
                    if not is_valid:
                        st.error(msg)
                        st.stop()
                        
                    # Merge with PRS results by Sample_ID
                    merged = pd.merge(prs_df, pheno_df, on='Sample_ID', how='inner')
                    
                    n_prs = len(prs_df)
                    n_matched = len(merged)
                    
                    if n_matched > 0:
                        phenotype = merged['PHENO'].values
                        prs_df = merged.drop(columns=['PHENO'])
                        
                        if n_matched < n_prs:
                            st.warning(f"Only {n_matched} of {n_prs} PRS samples have phenotype values. ML will use {n_matched} matched samples.")
                        else:
                            st.info(f"Merged {n_matched} samples with phenotype data.")
                    else:
                        st.error("No matching Sample IDs between PRS results and phenotype file. Check that IDs match.")
                        st.stop()
                        
                elif use_demo_mode:
                    phenotype = generate_mock_phenotype(n_samples=len(prs_df), binary=ml_is_binary)
                    st.warning("**Demo Mode**: Using randomly generated phenotype. Results are for testing only.")
                else:
                    st.error("Please upload a Phenotype File or enable Demo Mode to proceed.")
                    st.stop()
                
                # --- Train ML models ---
                with st.spinner("Training ML Models..."):
                    X = prs_df[selected_features]
                    y = phenotype
                    sample_ids = prs_df['Sample_ID'].values if 'Sample_ID' in prs_df.columns else None
                    
                    ml_df, predictions, y_test, test_indices = train_ml_models(X, y, selected_ml, is_binary=ml_is_binary, sample_ids=sample_ids)
                    
                    st.session_state.ml_results = ml_df
                    st.session_state.phenotype = phenotype
                    st.session_state.ml_prs_df = prs_df
                    st.session_state.ml_predictions = predictions
                    st.session_state.ml_y_test = y_test
                    st.session_state.ml_test_indices = test_indices
                    
                    # Auto-save ML results to disk
                    import os
                    results_dir = "/app/results" if os.path.exists("/app") else "./results"
                    os.makedirs(results_dir, exist_ok=True)
                    ml_save_path = os.path.join(results_dir, "ml_results.csv")
                    ml_df.to_csv(ml_save_path)
                    
                st.success(f"ML Training Complete. Results saved to `{ml_save_path}`")
                st.rerun()

    # ================================================================
    # TABS: Results Display
    # ================================================================
    st.divider()
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "Single Ancestry", 
        "Multi Ancestry", 
        "ML Predictors", 
        "Comparison & Visualization", 
        "Results & Export"
    ])
    
    with tab1:
        st.header("Single Ancestry Analysis")
        if config["mode"] == "Single Ancestry":
            if st.session_state.prs_results is not None:
                st.write(f"**Target Ancestry:** {config['ancestry']}")
                st.dataframe(st.session_state.prs_results.head(20))
            else:
                st.info("Run the PRS Calculation to view Single Ancestry results.")
        else:
            st.info("Currently in Multi Ancestry mode. Switch mode in sidebar to use Single Ancestry features.")

    with tab2:
        st.header("Multi Ancestry Integration")
        if config["mode"] == "Multi Ancestry":
            if st.session_state.prs_results is not None:
                st.write(f"**Target Ancestries Integrated:** {', '.join(config['ancestry']) if isinstance(config['ancestry'], list) else config['ancestry']}")
                st.dataframe(st.session_state.prs_results.head(20))
            else:
                st.info("Run the PRS Calculation to view Multi Ancestry results.")
        else:
            st.info("Currently in Single Ancestry mode. Switch mode in sidebar to use Multi Ancestry features.")

    with tab3:
        st.header("Machine Learning Predictors")
        if st.session_state.ml_results is not None:
            st.subheader("Model Performance")
            st.dataframe(st.session_state.ml_results)
            st.markdown("Models are trained using the PRS score columns as features.")
        elif st.session_state.prs_results is not None:
            st.info("PRS results are ready. Upload a phenotype file and click 'Run ML Prediction' above to train ML models.")
        else:
            st.info("Please run PRS Calculation first, then use 'Run ML Prediction' to train models.")

    with tab4:
        st.header("Comparison & Visualization")
        if st.session_state.prs_results is not None:
            prs_cols = [c for c in st.session_state.prs_results.columns if c != 'Sample_ID']
            
            col_a, col_b = st.columns(2)
            with col_a:
                st.subheader("PRS Distributions")
                hist_data = [st.session_state.prs_results[c].dropna() for c in prs_cols]
                fig_hist = ff.create_distplot(hist_data, prs_cols, show_hist=False)
                fig_hist.update_layout(margin=dict(l=0, r=0, t=30, b=0))
                st.plotly_chart(fig_hist, use_container_width=True)
                
            with col_b:
                st.subheader("Correlation Heatmap")
                corr = st.session_state.prs_results[prs_cols].corr()
                fig_corr = px.imshow(corr, text_auto=True, color_continuous_scale='RdBu_r', aspect="auto")
                fig_corr.update_layout(margin=dict(l=0, r=0, t=30, b=0))
                st.plotly_chart(fig_corr, use_container_width=True)
                
            st.subheader("PRS vs Phenotype")
            if st.session_state.phenotype is not None and st.session_state.ml_prs_df is not None and len(prs_cols) > 0:
                selected_prs_plot = st.selectbox("Select PRS Method for Scatter Plot", prs_cols)
                plot_df = pd.DataFrame({
                    'PRS Score': st.session_state.ml_prs_df[selected_prs_plot].values,
                    'Phenotype': st.session_state.phenotype
                })
                fig_scatter = px.scatter(plot_df, x='PRS Score', y='Phenotype', trendline="ols", opacity=0.6)
                st.plotly_chart(fig_scatter, use_container_width=True)
            else:
                st.info("Upload phenotype and run ML Prediction to see PRS vs Phenotype scatter plots.")
                
            if hasattr(st.session_state, 'ml_predictions') and st.session_state.ml_predictions is not None:
                st.subheader("ML Predictions vs True Phenotype (Test Set)")
                selected_model = st.selectbox("Select ML Model", list(st.session_state.ml_predictions.keys()))
                
                y_true = st.session_state.ml_y_test
                y_pred = st.session_state.ml_predictions[selected_model]
                
                plot_ml_df = pd.DataFrame({
                    'True Phenotype': y_true,
                    'Predicted Value': y_pred
                })
                
                fig_ml = px.scatter(plot_ml_df, x='True Phenotype', y='Predicted Value', trendline="ols", opacity=0.6)
                st.plotly_chart(fig_ml, use_container_width=True)

        else:
            st.info("Run the PRS Calculation to generate interactive visualizations.")

    with tab5:
        st.header("Results & Export")
        if st.session_state.prs_results is not None:
            st.subheader("PRS Scores")
            st.dataframe(st.session_state.prs_results)
            
            st.subheader("Download Options")
            col_c, col_d = st.columns(2)
            with col_c:
                csv_data = export_to_csv(st.session_state.prs_results)
                st.download_button(
                    label="Download PRS as CSV",
                    data=csv_data,
                    file_name="prs_results.csv",
                    mime="text/csv",
                    use_container_width=True
                )
            with col_d:
                excel_data = export_to_excel(st.session_state.prs_results)
                st.download_button(
                    label="Download PRS as Excel",
                    data=excel_data,
                    file_name="prs_results.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    use_container_width=True
                )
            
            # ML Results export
            if st.session_state.ml_results is not None:
                st.divider()
                st.subheader("ML Model Performance")
                st.dataframe(st.session_state.ml_results)
                ml_csv = export_to_csv(st.session_state.ml_results)
                st.download_button(
                    label="Download ML Results as CSV",
                    data=ml_csv,
                    file_name="ml_results.csv",
                    mime="text/csv",
                    use_container_width=True
                )
        else:
            st.info("No results to export. Run the PRS Calculation first.")

if __name__ == "__main__":
    main()
