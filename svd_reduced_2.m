reg_log;

% Filtra a matriz reduzida Ar para conter apenas os genes alvo fornecidos
fid = fopen('genes_nao_sexuais.txt', 'r');
if fid == -1
    error('Não foi possível abrir o arquivo genes_nao_sexuais.txt. Execute filter_genes.py primeiro.');
end
target_genes_raw = textscan(fid, '%s');
fclose(fid);
target_genes = target_genes_raw{1};

% Obtendo os nomes dos marcadores presentes em Ar (iSel vem de reg_log.m)
selected_genes = geneNames(iSel);

% Criando uma máscara binária para saber quais das linhas pertencem ao nosso subset
mask = ismember(selected_genes, target_genes);

% Aplica o filtro em Ar
Ar = Ar(mask, :);
disp('Nova dimensão de Ar após filtragem:');
disp(size(Ar));

% SVD - Matriz ultra reduzida
[U, S, V] = svd(Ar);
s = diag(S);
disp('Singular values:');
disp(s);
s = s*(1/sum(s));

figure;
plot(s);
hold on;
set(gcf, 'color', 'none');
plot(s, '*');
title('Valores singulares normalizados');
print(gcf, '-dpng', 'images/singular_values_reduced_2.png', '-r300')
clf

Aux2 = S*V';

x = Aux2(1, :);       % uso do padrão 1
y = Aux2(2, :);       % uso do padrão 2 
z = Aux2(3, :);       % uso do padrão 3

% Plotando classe 0 em verde e classe 1 em vermelho
figure;
hold on;
idx0 = find(b == 0);
idx1 = find(b == 1);

plot3(x(idx0), y(idx0), z(idx0), 'og', 'MarkerFaceColor', 'g', 'MarkerSize', 4, 'DisplayName', 'Class 0 (Healthy)');
plot3(x(idx1), y(idx1), z(idx1), 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 4, 'DisplayName', 'Class 1 (FMA)');

set(gcf, 'color', 'none');
set(gcf, 'Position', [100, 100, 800, 600]);

view(45, 30);  % Set viewing angle to show all 3 axes
title('SVD - Primeiros 3 componentes principais');
xlabel('Padrão 1');
ylabel('Padrão 2');
zlabel('Padrão 3');
legend('Location', 'best');

% Save a still image of the 3D SVD plot before creating the GIF
print(gcf, '-dpng', '-r300', 'images/svd_results_reduced_2.png')
filename = 'images/svd_results_reduced_2.gif';

% Force the axes to stay constant and prevent bouncing during rotation
axis tight manual;
axis vis3d;

for k = 1:72
    view(k*5, 30);   % rotate around the plot
    drawnow;

    frame = getframe(gcf);
    im = frame2im(frame);
    [A_img, map] = rgb2ind(im, 256);

    if k == 1
        imwrite(A_img, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.08);
    else
        imwrite(A_img, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.08);
    end
end

hold off;