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
    page_title="PRS Dashboard - Multi & Single Ancestry + ML Predictors",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- Main Application Logic ---
def main():
    st.title("🧬 PRS Dashboard - Multi & Single Ancestry + ML Predictors")
    st.markdown("A complete platform for polygenic risk score calculation, multi-ancestry integration, and machine learning prediction.")
    
    # Render Sidebar and get configs
    config = render_sidebar()
    
    # Initialize Session State
    if 'prs_results' not in st.session_state:
        st.session_state.prs_results = None
    if 'ml_results' not in st.session_state:
        st.session_state.ml_results = None
    if 'phenotype' not in st.session_state:
        st.session_state.phenotype = None

    # ================================================================
    # SECTION 1: PRS Calculation
    # ================================================================
    st.header("📊 PRS Calculation")
    
    col_prs_btn, col_prs_status = st.columns([1, 3])
    with col_prs_btn:
        run_prs_btn = st.button("🧬 Run PRS Calculation", type="primary", use_container_width=True)
    with col_prs_status:
        if st.session_state.prs_results is not None:
            n_samples = len(st.session_state.prs_results)
            n_scores = len([c for c in st.session_state.prs_results.columns if c != 'Sample_ID'])
            st.success(f"✅ PRS results available: {n_samples} samples, {n_scores} score column(s)")
        else:
            st.info("No PRS results yet. Upload data and click 'Run PRS Calculation'.")

    if run_prs_btn:
        if not config["selected_methods"]:
            st.warning("⚠️ Please select at least one PRS method from the sidebar.")
        else:
            with st.spinner("Processing Genotypes and Calculating PRS..."):
                import os
                
                data_dir = "/app/data" if os.path.exists("/app") else "./data"
                os.makedirs(data_dir, exist_ok=True)
                
                # --- Validate and save GWAS files ---
                gwas_paths = []
                gwas_pops = []
                gwas_ns = []
                if config["gwas_info"] and len(config["gwas_info"]) > 0:
                    for info in config["gwas_info"]:
                        f = info["file"]
                        path = os.path.join(data_dir, f.name)
                        with open(path, "wb") as out_f:
                            out_f.write(f.getbuffer())
                        gwas_paths.append(path)
                        gwas_pops.append(info["pop"])
                        gwas_ns.append(str(info["n_gwas"]))
                else:
                    st.error("❌ Please upload at least one GWAS summary statistics file.")
                    st.stop()
                
                # --- Validate and save target genotype ---
                target_prefix = "mock_target"
                if config["target_data"] is not None and config["target_data"] != []:
                    exts = [f.name.split('.')[-1] for f in config["target_data"]]
                    if not all(ext in exts for ext in ['bed', 'bim', 'fam']):
                        st.error("❌ Please upload all three PLINK binary files: .bed, .bim, .fam")
                        st.stop()
                    
                    for f in config["target_data"]:
                        path = os.path.join(data_dir, f.name)
                        with open(path, "wb") as out_f:
                            out_f.write(f.getbuffer())
                    
                    bed_file = [f.name for f in config["target_data"] if f.name.endswith('.bed')][0]
                    target_prefix = os.path.join(data_dir, bed_file.rsplit('.bed', 1)[0])
                else:
                    st.error("❌ Please upload Target Genotype Data (.bed/.bim/.fam).")
                    st.stop()
                
                # --- Validation data (optional) ---
                val_prefix = None
                if config["val_data"] is not None and config["val_data"] != []:
                    exts = [f.name.split('.')[-1] for f in config["val_data"]]
                    if not all(ext in exts for ext in ['bed', 'bim', 'fam']):
                        st.error("❌ Please upload all three PLINK binary files for validation: .bed, .bim, .fam")
                        st.stop()
                    for f in config["val_data"]:
                        path = os.path.join(data_dir, "val_" + f.name)
                        with open(path, "wb") as out_f:
                            out_f.write(f.getbuffer())
                    bed_file = [f.name for f in config["val_data"] if f.name.endswith('.bed')][0]
                    val_prefix = os.path.join(data_dir, "val_" + bed_file.rsplit('.bed', 1)[0])
                    
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
                    st.error("🚨 No PRS results were generated. Please check the tool errors above.")
                    st.stop()
                
                if 'Sample_ID' not in prs_df.columns:
                    prs_df.insert(0, 'Sample_ID', [f"IID_{i}" for i in range(len(prs_df))])
                
                # Store PRS results
                st.session_state.prs_results = prs_df
                # Reset ML results since PRS data changed
                st.session_state.ml_results = None
                st.session_state.phenotype = None
                
            st.success(f"✅ PRS Calculation Complete! {len(prs_df)} samples scored.")
            st.rerun()

    # ================================================================
    # SECTION 2: ML Prediction
    # ================================================================
    st.divider()
    st.header("🤖 ML Prediction")
    
    if st.session_state.prs_results is None:
        st.warning("⚠️ Please run PRS Calculation first before training ML models.")
    else:
        # Phenotype upload for ML (in main area, not sidebar)
        ml_col1, ml_col2 = st.columns(2)
        with ml_col1:
            ml_pheno_file = st.file_uploader(
                "📄 Upload Phenotype File for ML Training",
                help="File with columns: FID IID Phenotype (or IID Phenotype). Required for real ML evaluation.",
                key="ml_pheno_upload"
            )
        with ml_col2:
            ml_is_binary = st.checkbox(
                "Is phenotype binary (case/control)?",
                value=config["is_binary"],
                key="ml_binary_check"
            )
            use_demo_mode = st.checkbox(
                "🧪 Demo mode (use random phenotype)",
                value=False,
                help="If enabled, generates a random phenotype for testing purposes only. Results will NOT be scientifically meaningful."
            )
        
        # Show current PRS columns available
        score_cols = [c for c in st.session_state.prs_results.columns if c != 'Sample_ID']
        st.caption(f"Available PRS features for ML: `{', '.join(score_cols)}`")
        
        run_ml_btn = st.button("🚀 Run ML Prediction", type="primary", use_container_width=False)
        
        if run_ml_btn:
            prs_df = st.session_state.prs_results.copy()
            phenotype = None
            
            # --- Load phenotype ---
            if ml_pheno_file is not None:
                try:
                    pheno_df = pd.read_csv(ml_pheno_file, sep='\s+|,', engine='python')
                    # Phenotype file: FID IID PHENO or just IID PHENO
                    if len(pheno_df.columns) >= 3:
                        pheno_df = pheno_df.iloc[:, [1, -1]]  # IID + last col
                    else:
                        pheno_df = pheno_df.iloc[:, [0, -1]]  # IID + last col
                    pheno_df.columns = ['Sample_ID', 'PHENO']
                    
                    # Merge with PRS results
                    merged = pd.merge(prs_df, pheno_df, on='Sample_ID', how='inner')
                    if len(merged) > 0:
                        phenotype = merged['PHENO'].values
                        prs_df = merged.drop(columns=['PHENO'])
                        st.info(f"📊 Merged {len(merged)} samples with phenotype data.")
                    else:
                        st.error("❌ No matching Sample IDs between PRS results and phenotype file!")
                        st.stop()
                except Exception as e:
                    st.error(f"❌ Error reading phenotype file: {e}")
                    st.stop()
                    
            elif use_demo_mode:
                phenotype = generate_mock_phenotype(n_samples=len(prs_df), binary=ml_is_binary)
                st.warning("🧪 **Demo Mode**: Using randomly generated phenotype. Results are for testing only!")
            else:
                st.error("❌ Please upload a Phenotype File or enable Demo Mode to proceed.")
                st.stop()
            
            # --- Train ML models ---
            with st.spinner("Training ML Models..."):
                X = prs_df[score_cols]
                y = phenotype
                
                selected_ml = config["selected_ml"]
                if not selected_ml:
                    st.warning("⚠️ No ML models selected. Please select at least one model from the sidebar.")
                    st.stop()
                
                ml_df, predictions = train_ml_models(X, y, selected_ml, is_binary=ml_is_binary)
                
                st.session_state.ml_results = ml_df
                st.session_state.phenotype = phenotype
                
            st.success("✅ ML Training Complete!")
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
            st.info("🤖 PRS results are ready. Upload a phenotype file and click 'Run ML Prediction' above to train ML models.")
        else:
            st.info("📋 Please run PRS Calculation first, then use 'Run ML Prediction' to train models.")

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
                
            st.subheader("PRS vs Phenotype (Scatter Plot)")
            if st.session_state.phenotype is not None and len(prs_cols) > 0:
                selected_prs_plot = st.selectbox("Select PRS Method for Scatter Plot", prs_cols)
                plot_df = pd.DataFrame({
                    'PRS Score': st.session_state.prs_results[selected_prs_plot],
                    'Phenotype': st.session_state.phenotype[:len(st.session_state.prs_results)]
                })
                fig_scatter = px.scatter(plot_df, x='PRS Score', y='Phenotype', trendline="ols", opacity=0.6)
                st.plotly_chart(fig_scatter, use_container_width=True)
            else:
                st.info("Upload phenotype and run ML Prediction to see PRS vs Phenotype scatter plots.")

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
                    label="📥 Download PRS as CSV",
                    data=csv_data,
                    file_name="prs_results.csv",
                    mime="text/csv",
                    use_container_width=True
                )
            with col_d:
                excel_data = export_to_excel(st.session_state.prs_results)
                st.download_button(
                    label="📥 Download PRS as Excel",
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
                    label="📥 Download ML Results as CSV",
                    data=ml_csv,
                    file_name="ml_results.csv",
                    mime="text/csv",
                    use_container_width=True
                )
        else:
            st.info("No results to export. Run the PRS Calculation first.")

if __name__ == "__main__":
    main()
