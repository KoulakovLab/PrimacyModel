close all
clear
clc

%% Load and plot the data

[corr_mx, corr_mz, DX, DZ] = load_datasets('no_shuffle');

%% Correlations

figure, imagesc(triu(corr_mx, 1)), axis square, caxis([-0.2, 0.2]);
colormap(white_colormap), colorbar
title('FlyEM dataset')
xlabel('Olfactory receptors')
ylabel('Olfactory receptors')

figure, imagesc(triu(corr_mz, 1)), axis square, caxis([-0.2, 0.2]);
colormap(white_colormap), colorbar
title('FAFB dataset')
xlabel('Olfactory receptors')
ylabel('Olfactory receptors')

%% Correlations of correlations

el_corr_mx = triu(corr_mx, 1);
el_corr_mx = el_corr_mx(el_corr_mx(:) ~= 0);
el_corr_mz = triu(corr_mz, 1);
el_corr_mz = el_corr_mz(el_corr_mz(:) ~= 0);

c = corr(el_corr_mx(:), el_corr_mz(:));

figure, scatter(el_corr_mx(:), el_corr_mz(:), '.')
axis([-.2 .2 -.2 .2]), axis square, box on
title(sprintf('Correlation of projection patterns %.2f', c))
xlabel('FlyEM dataset')
ylabel('FAFB dataset')

%shuffled
[~, corr_mz_shuffled, ~, ~] = load_datasets('shuffle');

el_corr_mz = triu(corr_mz_shuffled, 1);
el_corr_mz = el_corr_mz(el_corr_mz(:) ~= 0);

c_shuffled = corr(el_corr_mx(:), el_corr_mz(:));

%bootstrap (disabled here for speedup)
el_corr_mz_bs = triu(corr_mz, 1);
el_corr_mz_bs = el_corr_mz_bs(el_corr_mz_bs(:) ~= 0);

c_bootstrap = corr(el_corr_mx(:), el_corr_mz_bs(:));

figure, plot([1, 2], [c_bootstrap, c_shuffled], 'o')
set(findobj(gca,'type','line'),'linew',3)
axis([0, 3, -0.1 0.6]), axis square
xticks([1,2])
xticklabels({'Bootstrap', 'Shuffle'})
ylabel('Pearson correlation')
