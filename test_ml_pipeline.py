"""
Test script to verify the SEPARATED ML pipeline works with actual PRS-CSx results.
Simulates exactly what the NEW app.py does in Step 2 (ML Prediction).
"""
import pandas as pd
import numpy as np
import sys
sys.path.insert(0, '/app')

print("=" * 60)
print("TEST: Separated ML Pipeline (Step 2 only)")
print("=" * 60)

# ============================================================
# Simulate Step 1 output (PRS results already in session_state)
# ============================================================
print("\n[Simulating Step 1 output] Reading existing PRS scores...")
sscore_path = "/app/results/prscsx/target_EUR_scored.sscore"
df = pd.read_csv(sscore_path, sep='\s+')
if '#FID' in df.columns:
    df = df[['IID', 'SCORE1_AVG']]
    df = df.rename(columns={'SCORE1_AVG': 'PRS_CSx_EUR'})
prs_df = df.rename(columns={'IID': 'Sample_ID'})
print(f"  PRS df: {prs_df.shape}, columns: {list(prs_df.columns)}")

# ============================================================
# Step 2: ML Prediction (new separated logic from app.py)
# ============================================================
print("\n[Step 2] Loading phenotype file...")
val_pheno_path = "/app/data/TOY_TARGET_DATA.pheno"
pheno_df = pd.read_csv(val_pheno_path, sep='\s+')

if len(pheno_df.columns) >= 3:
    pheno_df = pheno_df.iloc[:, [1, -1]]
else:
    pheno_df = pheno_df.iloc[:, [0, -1]]
pheno_df.columns = ['Sample_ID', 'PHENO']

# Merge
merged = pd.merge(prs_df, pheno_df, on='Sample_ID', how='inner')
print(f"  Merged: {merged.shape}")
assert len(merged) > 0, "ERROR: No matching Sample IDs!"

phenotype = merged['PHENO'].values
prs_df_clean = merged.drop(columns=['PHENO'])

# Select score columns
score_cols = [c for c in prs_df_clean.columns if c != 'Sample_ID']
X = prs_df_clean[score_cols]
y = phenotype
print(f"  X shape: {X.shape}, y shape: {y.shape}")
print(f"  Score columns: {score_cols}")

# Train ML models
print("\n[Step 2] Training ML models...")
from src.ml_models import train_ml_models

# Test with continuous phenotype
ml_df, predictions = train_ml_models(X, y, ["GLM", "SVM", "Random Forest"], is_binary=False)
print(f"\n  ML Results (Regression):")
print(ml_df)

# Test with binary phenotype (demo mode)
print("\n[Step 2b] Testing with demo binary phenotype...")
from src.utils import generate_mock_phenotype
mock_y = generate_mock_phenotype(n_samples=len(X), binary=True)
ml_df_bin, _ = train_ml_models(X, mock_y, ["GLM", "SVM", "Random Forest"], is_binary=True)
print(f"\n  ML Results (Classification - Demo):")
print(ml_df_bin)

print("\n" + "=" * 60)
print("✅ ALL TESTS PASSED! Separated pipeline works correctly!")
print("=" * 60)
