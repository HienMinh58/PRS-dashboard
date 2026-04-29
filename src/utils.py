import pandas as pd
import numpy as np
from sklearn.metrics import r2_score, roc_auc_score, mean_squared_error
from scipy.stats import pearsonr
import io

def generate_mock_phenotype(n_samples=500, binary=False):
    """Generates mock phenotypes for testing."""
    np.random.seed(42)
    if binary:
        return np.random.binomial(1, 0.3, n_samples)
    else:
        return np.random.normal(0, 1, n_samples)

def calculate_metrics(y_true, y_pred, is_binary=False):
    """Calculate performance metrics."""
    metrics = {}
    try:
        if is_binary:
            metrics['AUC'] = roc_auc_score(y_true, y_pred)
        else:
            metrics['R2'] = r2_score(y_true, y_pred)
            metrics['MSE'] = mean_squared_error(y_true, y_pred)
            corr, pval = pearsonr(y_true, y_pred)
            metrics['Correlation'] = corr
            metrics['P-value'] = pval
    except Exception as e:
        metrics['Error'] = str(e)
    return metrics

def export_to_csv(df):
    """Convert dataframe to CSV bytes for download."""
    return df.to_csv(index=False).encode('utf-8')

def export_to_excel(df):
    """Convert dataframe to Excel bytes for download."""
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        df.to_excel(writer, index=False, sheet_name='Results')
    return output.getvalue()
