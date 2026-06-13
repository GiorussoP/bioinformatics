import pandas as pd
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, accuracy_score
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.base import BaseEstimator, TransformerMixin
import warnings
warnings.filterwarnings('ignore')

# Classe customizada que encapsula o algoritmo 'resolve' do MATLAB sem Data Leakage
class RidgeSystemSelector(BaseEstimator, TransformerMixin):
    def __init__(self, n_features_to_select=33):
        self.n_features_to_select = n_features_to_select
        self.selected_indices_ = None
        
    def fit(self, X, y):
        # X no sklearn possui dimensões (amostras, genes)
        m, n = X.shape
        Im = sp.eye(m, format='csr')
        In = sp.eye(n, format='csr')
        X_sparse = sp.csr_matrix(X)
        
        # Reconstrói a matriz aumentada exata do MATLAB M = [Im, -A; A', In]
        M = sp.bmat([[Im, -X_sparse], [X_sparse.T, In]], format='csr')
        
        nb = np.zeros(m + n)
        nb[:m] = -y
        
        # Resolve o sistema esparso linear
        x_sol = spsolve(M, nb)
        w = x_sol[m:] # Pesos calculados para cada gene
        
        # Seleciona as k melhores características baseadas na magnitude de w
        self.selected_indices_ = np.argsort(np.abs(w))[::-1][:self.n_features_to_select]
        return self
        
    def transform(self, X):
        return X[:, self.selected_indices_]

# 1. Carregar marcadores autossômicos válidos
with open('genes_nao_sexuais.txt', 'r') as f:
    genes_validos = [line.strip() for line in f if line.strip()]

# 2. Carregar a base completa gerada pelo MATLAB
print("Carregando a matriz genômica completa...")
df = pd.read_csv('expressao_genes_completa.csv')

# 3. Aplicar filtro cromossômico global seguro (não causa leakage)
df_filtrado = df[df['Marcador'].astype(str).isin(genes_validos)].set_index('Marcador')
X = df_filtrado.T.values
pacientes = df_filtrado.columns

# 4. Criar os rótulos dinamicamente baseados nos nomes das colunas
y = np.array([1 if 'FMA' in pac else 0 for pac in pacientes])

print(f"Dimensões prontas para o modelo: {X.shape[0]} amostras e {X.shape[1]} genes.")

# 5. Divisão Holout (Treino/Teste) rigorosa
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42, stratify=y)

# 6. Construção do Pipeline de Machine Learning Robusto
n_components = 5
pipeline = Pipeline([
    ('scaler', StandardScaler()),                                # Normalização por fold
    ('selector', RidgeSystemSelector(n_features_to_select=33)),  # Seleção matemática por fold (SEM LEAKAGE!)
    ('svd', TruncatedSVD(n_components=n_components, random_state=42)),
    ('clf', LogisticRegression(class_weight='balanced', random_state=42))
])

# Ajustar e avaliar no Hold-out
pipeline.fit(X_train, y_train)
y_pred_svd = pipeline.predict(X_test)
acc_svd = accuracy_score(y_test, y_pred_svd)

print("\n--- Resultados da Classificação Segura (Hold-out Test) ---")
print(f"Acurácia no conjunto de teste: {acc_svd * 100:.2f}%\n")
print(classification_report(y_test, y_pred_svd, target_names=['Healthy', 'FMA']))

# 7. Validação Cruzada de 5 Folds com blindagem total contra vazamento de dados
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
scores = cross_val_score(pipeline, X, y, cv=cv, scoring='accuracy')

print("--- Resultados da Validação Cruzada (5-Fold CV Robusto) ---")
print(f"Acurácias por fold: {[f'{s*100:.2f}%' for s in scores]}")
print(f"Acurácia Média Final: {scores.mean() * 100:.2f}% (+/- {scores.std() * 100:.2f}%)")