% export_completo.m
clear; clc;
load_data; % Carrega A, b e geneNames originais do .mat

[n_genes, n_patients] = size(A);

% Criar os cabeçalhos das colunas identificando as classes dos pacientes
header = cell(1, n_patients + 1);
header{1} = 'Marcador';

for j = 1:n_patients
    if b(j) == 1
        header{j+1} = sprintf('PacienteFMA_%d', j);
    else
        header{j+1} = sprintf('PacienteHealthy_%d', j);
    end
end

fprintf('Convertendo a matriz completa A (%d genes x %d pacientes)...\n', n_genes, n_patients);

% Transforma em tabela para exportação direta e rápida
T = array2table(A);
T.Properties.VariableNames = header(2:end);

% Adiciona a coluna com os símbolos dos genes no início
Marcador = geneNames;
T = [table(Marcador) T];

fprintf('Escrevendo expressao_genes_completa.csv (isso pode levar alguns segundos)...\n');
writetable(T, 'expressao_genes_completa.csv');
fprintf('Pronto! Matriz completa exportada com sucesso.\n');