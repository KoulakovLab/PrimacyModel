close all
clear
clc

%% Load the data

[~, ~, DX, DZ] = load_datasets('no_shuffle');
DX_data{1} = DX;
DX_data{2} = DZ;
[~, ~, DX, DZ] = load_datasets('shuffle');
DX_data{3} = DX;
DX_data{4} = DZ;

%% PCA
cmap = lines;
    
for i = 1 : 2
    DX = DX_data{i};
    [~, explained_pca] = cmdscale(DX);
    explained_pca = explained_pca(1 : end - 1);
    explained_pca = explained_pca / sum(explained_pca) * 100;
    figure(1), semilogx(sort(explained_pca, 'descend'), '-o', ...
        'linewidth', 3, 'color', cmap(i, :), 'markerfacecolor', cmap(i, :)); hold on
    figure(2), semilogx(cumsum(sort(explained_pca, 'descend')), '-o', ...
        'linewidth', 3, 'color', cmap(i, :), 'markerfacecolor', cmap(i, :)); hold on
end

for i = 3 : 4
    DX = DX_data{i};
    [~, explained_pca] = cmdscale(DX);
    explained_pca = explained_pca(1 : end - 1);
    explained_pca = explained_pca / sum(explained_pca) * 100;
    figure(1), semilogx(sort(explained_pca, 'descend'), ':o', ...
        'linewidth', 3, 'color', cmap(i-2, :), 'markerfacecolor', cmap(i-2, :))
    figure(2), semilogx(cumsum(sort(explained_pca, 'descend')), ':o', ...
        'linewidth', 3, 'color', cmap(i-2, :), 'markerfacecolor', cmap(i-2, :))
end

figure(1), xlabel('Dimension #'), ylabel('Variance explained'), title('PCA')
legend("FlyEM", "FAFB", "FlyEM shfld", "FAFB shfld")
figure(2), xlabel('Dimension #'), ylabel('Cumulative variance explained')

%% PCA - LowD structure

[score1, ~] = cmdscale(DX_data{1});
[score2, ~] = cmdscale(DX_data{2});

figure(5), scatter(score1(:, 1), score1(:, 2), 100, 'markeredgecolor', [0, 0.4470, 0.7410],  'markerfacecolor', [0, 0.4470, 0.7410], 'linewidth', 3), hold on
scatter(score2(:, 1), -score2(:, 2), 100, 'markeredgecolor', [0.8500, 0.3250, 0.0980], 'markerfacecolor', [0.8500, 0.3250, 0.0980], 'linewidth', 3)

for i = 1 : size(score1, 1)
    line([score1(i, 1), score2(i, 1)], [score1(i, 2), -score2(i, 2)], 'color', 'k', 'linewidth', 2)
end
legend('FlyEM', 'FAFB');
title('2D PCA embedding')
xlabel('PCA dimansion 1');
ylabel('PCA dimension 2');
box on

%% Isomap

options.verbose = 0;
options.display = 0;

DX_tmp = cell(4, 50);

for i = 1 : 4
    DX = DX_data{i};
    
    options.dims = size(DX, 1);

    explained_var_k = zeros(size(DX, 1) - 1, 10);

    for k = 2 : 50
        if i == 3
            [~, ~, DX, ~] = load_datasets('shuffle');
        elseif i == 4
            [~, ~, ~, DX] = load_datasets('shuffle');
        end
        [DX_isomap, residual_isomap, ~] = Isomap(DX, 'k', k, options);
        DX_isomap = DX_isomap.coords{1};
        [~, ~, ~, ~, explained_isomap, ~] = pca(DX_isomap');
        explained_var_k(:, k) = sort(explained_isomap, 'descend');
        
        DX_tmp{i, k} = DX_isomap;
    end
end

check_Isomap_residuals2;

%% Isomap - evaluation

k_max = 6;

for i = 1 : 4
    DX = DX_data{i};
    
    options.dims = size(DX, 1);

    if i == 3
        [~, ~, DX, ~] = load_datasets('shuffle');
    elseif i == 4
        [~, ~, ~, DX] = load_datasets('shuffle');
    end
    [DX_isomap, residual_isomap, ~] = Isomap(DX, 'k', k_max, options);
    DX_isomap = DX_isomap.coords{1};
    [~, ~, ~, ~, explained_isomap, ~] = pca(DX_isomap');
    explained_var_k = sort(explained_isomap, 'descend');

    %Plots (errorbars removed for speedup)
    if i < 3
        figure(3), semilogx(explained_var_k, '-o', ...
            'linewidth', 3, 'color', cmap(i, :), 'markerfacecolor', cmap(i, :)), hold on;
        figure(4), semilogx(cumsum(explained_var_k), '-o', ...
            'linewidth', 3, 'color', cmap(i, :), 'markerfacecolor', cmap(i, :)), hold on;
    else
        figure(3), semilogx(explained_var_k, ':o', ...
            'linewidth', 3, 'color', cmap(i-2, :), 'markerfacecolor', cmap(i-2, :));
        figure(4), semilogx(cumsum(explained_var_k), ':o', ...
            'linewidth', 3, 'color', cmap(i-2, :), 'markerfacecolor', cmap(i-2, :));
    end
    
    if i == 1
        DX_isomap1 = DX_isomap;
    elseif i == 2
        DX_isomap2 = DX_isomap;
    end

end

figure(3), xlabel('Dimension #'), ylabel('Variance explained')
figure(4), xlabel('Dimension #'), ylabel('Cumulative variance explained'), ...
    title('Isomap'), axis([0 100 0 100])

%% Isomap - LowD structure

score1 = DX_isomap1';
score2 = DX_isomap2';

figure(6), scatter(score1(:, 1), -score1(:, 2), 100, 'markeredgecolor', [0, 0.4470, 0.7410],  'markerfacecolor', [0, 0.4470, 0.7410], 'linewidth', 3), hold on
scatter(score2(:, 1), score2(:, 2), 100, 'markeredgecolor', [0.8500, 0.3250, 0.0980], 'markerfacecolor', [0.8500, 0.3250, 0.0980], 'linewidth', 3)

for i = 1 : size(score1, 1)
    line([score1(i, 1), score2(i, 1)], [-score1(i, 2), score2(i, 2)], 'color', 'k', 'linewidth', 2)
end
title('2D Isomap embedding')
xlabel('Isomap dimension 1');
ylabel('Isomap dimension 2');
box on
