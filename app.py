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

    # Top action button
    col1, col2 = st.columns([1, 4])
    with col1:
        run_btn = st.button("🚀 Run All Selected Methods", type="primary", use_container_width=True)
        
    if run_btn:
        if not config["selected_methods"]:
            st.warning("Please select at least one PRS method from the sidebar.")
        else:
            with st.spinner("Processing Genotypes and Calculating PRS..."):
                # Validate and Save Files to /app/data/
                import os
                
                # Check if data directory exists, if not use local temp for testing outside docker
                data_dir = "/app/data" if os.path.exists("/app") else "./data"
                os.makedirs(data_dir, exist_ok=True)
                
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
                    st.error("Please upload at least one GWAS summary statistics file.")
                    st.stop()
                
                target_prefix = "mock_target"
                if config["target_data"] is not None and config["target_data"] != []:
                    # Validate we have .bed, .bim, .fam
                    exts = [f.name.split('.')[-1] for f in config["target_data"]]
                    if not all(ext in exts for ext in ['bed', 'bim', 'fam']):
                        st.error("Please upload all three PLINK binary files: .bed, .bim, .fam")
                        st.stop()
                    
                    for f in config["target_data"]:
                        path = os.path.join(data_dir, f.name)
                        with open(path, "wb") as out_f:
                            out_f.write(f.getbuffer())
                    
                    # Prefix is the filename without extension of the .bed file
                    bed_file = [f.name for f in config["target_data"] if f.name.endswith('.bed')][0]
                    target_prefix = os.path.join(data_dir, bed_file.rsplit('.bed', 1)[0])
                else:
                    st.error("Please upload Target Genotype Data (.bed/.bim/.fam).")
                    st.stop()
                
                # Process validation data
                val_prefix = None
                if config["val_data"] is not None and config["val_data"] != []:
                    exts = [f.name.split('.')[-1] for f in config["val_data"]]
                    if not all(ext in exts for ext in ['bed', 'bim', 'fam']):
                        st.error("Please upload all three PLINK binary files for validation: .bed, .bim, .fam")
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

                # Run PRS Methods
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
                    st.error("🚨 Không có kết quả PRS nào được tạo ra. Vui lòng kiểm tra lỗi từ công cụ ở trên.")
                    st.stop()
                
                # For ML Predictors tab, we'll use the target PRS results (since validation was used to fit the combined PRS)
                # But ML predictors might need a phenotype for the target set (if training more models). 
                # If target phenotype is provided (in val_pheno or generated mock), we use it.
                if val_pheno_path:
                    try:
                        pheno_df = pd.read_csv(val_pheno_path, sep='\s+|,', engine='python')
                        st.session_state.phenotype = pheno_df.iloc[:, -1].values
                    except:
                        st.session_state.phenotype = generate_mock_phenotype(n_samples=len(prs_df), binary=config["is_binary"])
                else:
                    st.session_state.phenotype = generate_mock_phenotype(n_samples=len(prs_df), binary=config["is_binary"])
                
                if 'Sample_ID' not in prs_df.columns:
                    prs_df.insert(0, 'Sample_ID', [f"IID_{i}" for i in range(len(prs_df))])
                
                # Store results in session
                st.session_state.prs_results = prs_df

                
                # Run ML Models if selected
                if config["selected_ml"]:
                    st.toast("Training ML Models...")
                    X = prs_df[config["selected_methods"]]
                    y = st.session_state.phenotype
                    ml_df, _ = train_ml_models(X, y, config["selected_ml"], is_binary=config["is_binary"])
                    st.session_state.ml_results = ml_df
                    
            st.success("Analysis Complete!")

    # --- Tabs ---
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
                st.info("Run the pipeline to view Single Ancestry results.")
        else:
            st.info("Currently in Multi Ancestry mode. Switch mode in sidebar to use Single Ancestry features.")

    with tab2:
        st.header("Multi Ancestry Integration")
        if config["mode"] == "Multi Ancestry":
            if st.session_state.prs_results is not None:
                st.write(f"**Target Ancestries Integrated:** {', '.join(config['ancestry']) if isinstance(config['ancestry'], list) else config['ancestry']}")
                st.dataframe(st.session_state.prs_results.head(20))
            else:
                st.info("Run the pipeline to view Multi Ancestry results.")
        else:
            st.info("Currently in Single Ancestry mode. Switch mode in sidebar to use Multi Ancestry features.")

    with tab3:
        st.header("Machine Learning Predictors")
        if st.session_state.ml_results is not None:
            st.subheader("Model Performance")
            st.dataframe(st.session_state.ml_results)
            st.markdown("Models are trained using the selected PRS methods as features.")
        else:
            st.info("Select ML Predictors in the sidebar and run the pipeline to see performance metrics.")

    with tab4:
        st.header("Comparison & Visualization")
        if st.session_state.prs_results is not None:
            prs_cols = [c for c in st.session_state.prs_results.columns if c != 'Sample_ID']
            
            col_a, col_b = st.columns(2)
            with col_a:
                st.subheader("PRS Distributions")
                hist_data = [st.session_state.prs_results[c] for c in prs_cols]
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
                    'Phenotype': st.session_state.phenotype
                })
                fig_scatter = px.scatter(plot_df, x='PRS Score', y='Phenotype', trendline="ols", opacity=0.6)
                st.plotly_chart(fig_scatter, use_container_width=True)

        else:
            st.info("Run the pipeline to generate interactive visualizations.")

    with tab5:
        st.header("Results & Export")
        if st.session_state.prs_results is not None:
            st.dataframe(st.session_state.prs_results)
            
            st.subheader("Download Options")
            col_c, col_d = st.columns(2)
            with col_c:
                csv_data = export_to_csv(st.session_state.prs_results)
                st.download_button(
                    label="📥 Download as CSV",
                    data=csv_data,
                    file_name="prs_results.csv",
                    mime="text/csv",
                    use_container_width=True
                )
            with col_d:
                excel_data = export_to_excel(st.session_state.prs_results)
                st.download_button(
                    label="📥 Download as Excel",
                    data=excel_data,
                    file_name="prs_results.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    use_container_width=True
                )
        else:
            st.info("No results to export. Run the pipeline first.")

if __name__ == "__main__":
    main()
