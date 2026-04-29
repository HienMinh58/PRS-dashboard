import pandas as pd
from sklearn.svm import SVR, SVC
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, roc_auc_score

def train_ml_models(X, y, selected_models, is_binary=False):
    """
    Trains selected machine learning models using PRS scores as features.
    X: DataFrame of PRS scores
    y: Series/Array of phenotypes
    """
    results = {}
    predictions = {}
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    for model_name in selected_models:
        try:
            if model_name == "SVM":
                model = SVC(probability=True) if is_binary else SVR()
            elif model_name == "GLM":
                model = LogisticRegression() if is_binary else LinearRegression()
            elif model_name == "Random Forest":
                model = RandomForestClassifier(n_estimators=100, random_state=42) if is_binary else RandomForestRegressor(n_estimators=100, random_state=42)
            else:
                continue
                
            model.fit(X_train_scaled, y_train)
            
            if is_binary:
                y_pred_proba = model.predict_proba(X_test_scaled)[:, 1]
                score = roc_auc_score(y_test, y_pred_proba)
                metric_name = "AUC"
                predictions[model_name] = y_pred_proba
            else:
                y_pred = model.predict(X_test_scaled)
                score = r2_score(y_test, y_pred)
                metric_name = "R2"
                predictions[model_name] = y_pred
                
            results[model_name] = {
                "Metric": metric_name,
                "Score": round(score, 4),
                "Model Type": "Classification" if is_binary else "Regression"
            }
        except Exception as e:
            results[model_name] = {"Error": str(e)}
            
    return pd.DataFrame.from_dict(results, orient='index'), predictions
