import pandas as pd
import numpy as np
from sklearn.metrics import r2_score, roc_auc_score, mean_squared_error
from scipy.stats import pearsonr
import io

def generate_mock_phenotype(n_samples=500, binary=False):
    """
    Generates mock phenotype data for testing and demonstration.

    Args:
        n_samples (int): Number of samples to generate.
        binary (bool): If True, generates binary (0/1) traits. Otherwise, continuous.

    Returns:
        np.ndarray: Array of generated phenotype values.
    """
    np.random.seed(42)
    if binary:
        return np.random.binomial(1, 0.3, n_samples)
    else:
        return np.random.normal(0, 1, n_samples)

def calculate_metrics(y_true, y_pred, is_binary=False):
    """
    Calculates statistical performance metrics for PRS predictions.

    Args:
        y_true (array-like): Ground truth phenotype values.
        y_pred (array-like): Predicted PRS values or scores.
        is_binary (bool): Whether the phenotype is binary (affects metrics chosen).

    Returns:
        dict: Dictionary containing metrics like AUC, R2, MSE, Correlation, and P-value.
    """
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
    """
    Converts a pandas DataFrame to CSV bytes.

    Args:
        df (pd.DataFrame): The DataFrame to convert.

    Returns:
        bytes: CSV encoded as UTF-8 bytes, suitable for download.
    """
    return df.to_csv(index=False).encode('utf-8')

def export_to_excel(df):
    """
    Converts a pandas DataFrame to Excel (xlsx) bytes.

    Args:
        df (pd.DataFrame): The DataFrame to convert.

    Returns:
        bytes: Excel file content as bytes.
    """
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        df.to_excel(writer, index=False, sheet_name='Results')
    return output.getvalue()

