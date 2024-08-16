%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fig 5 Primacy 2024 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script will generate Fig. 5 F,G,H.
% The red box is drawn posthoc to bound the region of statistical significance defined by q < prm.Qval_thresh
% This can be obtained by running the script with prm.applyQmask = 1
% The figure without a significance mask is obtained with prm.applyQmask = 0

%%% 0
prm = struct();
prm.S_raw_impte = 'raw';        % raw vs imputed; imputing values in DoOR gives a stronger overlap result 
prm.m_topmost = [1 10 50];      % number of topmost overlaps coinsidered
prm.prange = 2:8;               % range of p numbers
prm.applyQmask = 0;             % view figure with q value mask
prm.Qval_thresh = 0.05;         % set q value threshold


%%% 1
[prm] = loadData(prm);


%%% 2
[anlyz] = analysis(prm);


%%% 3
plot(anlyz, prm);




%%% FNC

function [plt] = plot(anlyz, prm)
plt = struct();
f = figure;
for mi = 1 : length(prm.m_topmost)

    Y = anlyz.(strcat('M',num2str(prm.m_topmost(mi)))).Y;
    Yn = anlyz.(strcat('M',num2str(prm.m_topmost(mi)))).Yn;
    Q = anlyz.(strcat('M',num2str(prm.m_topmost(mi)))).Q;
    P = anlyz.(strcat('M',num2str(prm.m_topmost(mi)))).P;

    diffYYn = Y - mean(Yn, 3);
    
    subplot(1, length(prm.m_topmost), mi), 
    h = imagesc(diffYYn);
    xticks(0:2:max(prm.prange))
    xticklabels(0:2:max(prm.prange))
    yticks(2:2:max(prm.prange))
    yticklabels(2:2:max(prm.prange))
    xlabel('overlap')
    ylabel('p number')
    title(strcat('M', num2str(prm.m_topmost(mi))))
    set(gca, 'FontSize', 15, 'FontWeight', 'bold')
    axis square
    colormap(gca, redblueu)
    colorbar

    if prm.applyQmask
        mask = Q(:,1:max(prm.prange)+1) < prm.Qval_thresh;
        maskedData = diffYYn;
        maskedData(~mask) = NaN;
        alphaMask = double(mask); % Set masked values to opaque
        set(h, 'AlphaData', alphaMask); % Use AlphaData attribute for transparency
    end
end
end


function [anlyz] = analysis(prm)
anlyz = struct();
% compute overlap scores for each proportion of cells kept
for mi = 1 : length(prm.m_topmost)
    [anlyz] = overlapScore(anlyz, prm, prm.m_topmost(mi));
end
end


function [anlyz] = overlapScore(anlyz, prm, M)
minp = min(prm.prange);
maxp = max(prm.prange);
Y = zeros(maxp, maxp+1);
Y_null = zeros(maxp, maxp+1, prm.num_null);
P = Y;
Q = P;
W = P;

for pi = min(prm.prange) : maxp
    S = prm.S.(sprintf('p%d', pi)).S;
    X = S * prm.C';
    [X, num_zero2remove] = overlap_topM(X, M);
    f = overlap_spectrum_full(X);
    f(1) = f(1) - num_zero2remove;
    Y(pi, 1:length(f)) = f;
    
    for i1 = 1 : prm.num_null
        X_ = S * prm.Cn{i1}'; 
        [X_, num_zero2remove] = overlap_topM(X_, M);
        f = overlap_spectrum_full(X_);
        f(1) = f(1) - num_zero2remove;
        Y_null(pi, 1:length(f), i1) = f;
    end
    
    limits = find((Y(pi, 1:length(f))>0));
    W(pi, limits) = 1;
    
    [~, P(pi, limits)] = ...
        ttest2(squeeze(Y_null(pi, limits, :))', Y(pi, limits));

    Q(pi, limits) = mafdr(P(pi, limits), 'BHFDR', 1);
    
    if sum(isnan(P(pi, limits))) > 0            % For those NaN p
        P(pi, isnan(P(pi, :))) = 1;
        Q(pi, isnan(Q(pi, :))) = 1;
    end
    
    P(pi, setdiff(1:maxp+1, limits)) = 1;    
    Q(pi, setdiff(1:maxp+1, limits)) = 1;
end

anlyz.(strcat('M', num2str(M))).P = P;
anlyz.(strcat('M', num2str(M))).Q = Q;
anlyz.(strcat('M', num2str(M))).W = W;
anlyz.(strcat('M', num2str(M))).Y = Y;
anlyz.(strcat('M', num2str(M))).Yn = Y_null;
end


function [V, num_zero2remove] = overlap_topM(X, M)
[ii, ic] = sort(X', 'descend');
ii(M+1:end,:) = 0;
num_zero2remove = numel(ii(M+1:end));
V = zeros(size(X'));
for j = 1 : size(ii,2)
    V(ic(1:M,j), j) = ii(1:M, j);
end
V=V';
end


function [f] = overlap_spectrum_full(X)
overlaps = X(:);
f = zeros(max(overlaps)+1, 1);
for i = 1 : max(overlaps)+1
   f(i) = length(find(overlaps==i-1)); 
end
end


function [prm] = loadData(prm)
% load affinity
prm = loadAff(prm);
% load connectivity
prm = loadCon(prm);
end


function [prm] = loadCon(prm)
load connectivity_datasets_2021
prm.C = connectivity.FlyEM_HB.subset_37ors.kc_or_matrix;
prm.Cn = connectivity.FlyEM_HB.subset_37ors.null.kc_or_matrix;
prm.num_null = length(prm.Cn);
end


function [prm] = loadAff(prm)
load('DoOR_subset_used_primacyPLOS.mat')
prm.A = DoOR_subsets.subset1.affinity_data;
prm.A_lab = DoOR_subsets.subset1.labels(1,:);     
prm.S = DoOR_subsets.subset1.S_mats.(prm.S_raw_impte);
end
