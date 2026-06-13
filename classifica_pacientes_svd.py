import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, accuracy_score
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.pipeline import Pipeline
import warnings
warnings.filterwarnings('ignore')

# 1. Carregar os marcadores não-sexuais detectados automaticamente
with open('genes_nao_sexuais.txt', 'r') as f:
    marcadores_selecionados = [line.strip() for line in f if line.strip()]

# 2. Carregar a base de dados
print("Carregando os dados...")
df = pd.read_csv('expressao_genes_reduzidos.csv')

# 3. Filtrar e transpor os dados
df_filtrado = df[df['Marcador'].isin(marcadores_selecionados)].set_index('Marcador')
X = df_filtrado.T

# 4. Criar as labels
y = np.array([1 if 'FMA' in paciente else 0 for paciente in X.index])

# 5. Dividir em Treino e Teste
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42, stratify=y)

# 6. Escalonar
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# 7. Aplicar SVD (TruncatedSVD)
# Ao reduzir para um número ótimo de componentes (ex: 5), removemos o ruído ("colinearidade")
# que pode estar atrapalhando o modelo com as 13 features simultâneas.
n_components = 5
svd = TruncatedSVD(n_components=n_components, random_state=42)
X_train_svd = svd.fit_transform(X_train_scaled)
X_test_svd = svd.transform(X_test_scaled)

print(f"\nVariância explicada pelos {n_components} componentes da SVD: {svd.explained_variance_ratio_.sum()*100:.2f}%")

# 8. Regressão Logística baseada apenas nos componentes da SVD
# Usamos class_weight='balanced' para compensar qualquer leve desequilíbrio e C ajustado.
model_svd = LogisticRegression(class_weight='balanced', random_state=42)
model_svd.fit(X_train_svd, y_train)

# 9. Avaliação na partição de Teste
y_pred_svd = model_svd.predict(X_test_svd)
acc_svd = accuracy_score(y_test, y_pred_svd)

print("\n--- Resultados da Classificação com SVD (Hold-out Test) ---")
print(f"Número de componentes utilizados: {n_components}")
print(f"Acurácia no conjunto de teste: {acc_svd * 100:.2f}%\n")
print(classification_report(y_test, y_pred_svd, target_names=['Healthy', 'FMA']))

# 10. Validação Cruzada (K-Fold)
# Avaliação muito mais robusta para conjuntos pequenos
pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('svd', TruncatedSVD(n_components=n_components, random_state=42)),
    ('clf', LogisticRegression(class_weight='balanced', random_state=42))
])

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
scores = cross_val_score(pipeline, X, y, cv=cv, scoring='accuracy')

print("--- Resultados da Validação Cruzada (5-Fold CV) ---")
print(f"Acurácias por partição: {[f'{s*100:.2f}%' for s in scores]}")
print(f"Acurácia Média Robusta: {scores.mean() * 100:.2f}% (+/- {scores.std() * 100:.2f}%)")

