function [alpha, x] = resolve(A,b)
    [m,n] = size(A);
    Im = speye(m);
    In = speye(n);
    M = sparse(m+n,m+n);
    M = [Im,-A; A',In];
    nb = zeros(m+n,1);
    nb(1:m) = -b;
    x = M\nb;
    alpha = x(m+1:end);
end

load_data;

% Standarization of the original matrix A
A_mean = mean(A, 2);
A_sigma = std(A, 0, 2);
% Prevent division by zero
A_sigma(A_sigma == 0) = 0.00001;
A = (A - A_mean) ./ A_sigma;


w = resolve(A', b);
p = A' * w; % Raw linear predictions, roughly hovering around 0 and 1

% Plot the predictions
figure;
hold on;
idx0 = find(b == 0);
scatter(idx0, p(idx0), 'g', 'filled', 'DisplayName', 'Class 0 (Healthy)');
idx1 = find(b == 1);
scatter(idx1, p(idx1), 'r', 'filled', 'DisplayName', 'Class 1 (FMA)');
title('Predicted Probabilities');
xlabel('Patient Index');
ylabel('Probability (p)');
yline(0.5, '--k', 'Threshold = 0.5', 'HandleVisibility', 'off');
legend('Location', 'best');
hold off;
% Save the figure to a file
saveas(gcf, 'images/full_reg.png');


% This was the number with the best results (ironically)
n_attr = 42;
[~, pos] = sort(abs(w), 'descend');
iSel = pos(1:n_attr);
Ar = A(iSel, :);
wr = Ar' \ b;  % Fit weights using least squares, b is a column vector
pr = Ar' * wr; % Raw linear predictions

fprintf('Selected attributes (1-based indices):\n');
disp(iSel);
if exist('geneNames', 'var')
    selected_genes = geneNames(iSel);
    weights = w(iSel);
    T = table(iSel, selected_genes, weights, 'VariableNames', {'Index', 'GeneSymbol', 'Weight'});
    disp('Table of 30 Selected Attributes:');
    disp(T);
    
    % Export to CSV
    idx1_all = find(b == 1);
    idx0_all = find(b == 0);
    ordered_indices = [idx1_all; idx0_all];
    
    header = {'Marcador'};
    for i = 1:length(idx1_all)
        header{end+1} = sprintf('PacienteFMA%d', i);
    end
    for i = 1:length(idx0_all)
        header{end+1} = sprintf('PacienteHealthy%d', i);
    end
    
    Ar_ordered = Ar(:, ordered_indices);
    csv_data = cell(n_attr + 1, length(header));
    csv_data(1,:) = header;
    
    for i = 1:n_attr
        if iscell(selected_genes)
            csv_data{i+1, 1} = selected_genes{i};
        else
            csv_data{i+1, 1} = selected_genes(i);
        end
        for j = 1:length(ordered_indices)
            csv_data{i+1, j+1} = Ar_ordered(i, j);
        end
    end
    
    writecell(csv_data, 'expressao_genes_reduzidos.csv');
    fprintf('Exported reduced matrix with genes and patients to expressao_genes_reduzidos.csv\n');
end

% Plot the predictions for the reduced model
figure;
hold on;
idx0 = find(b == 0);
scatter(idx0, pr(idx0), 'g', 'filled', 'DisplayName', 'Class 0 (Healthy)');
idx1 = find(b == 1);
scatter(idx1, pr(idx1), 'r', 'filled', 'DisplayName', 'Class 1 (FMA)');
title('Predicted Probabilities - Attribute selection');
xlabel('Patient Index');
ylabel('Probability (p)');
yline(0.5, '--k', 'Threshold = 0.5', 'HandleVisibility', 'off');
legend('Location', 'best');
hold off;
% Save the figure to a file
saveas(gcf, 'images/reduced_reg.png');

% Save the reduced matrix
save('reduced_fibro_dataset.mat','Ar');
fprintf('Saved Ar to reduced_fibro_dataset.mat\n');
