c = zeros(1, 50) * NaN;

for i = 2 : 50
    A = DX_tmp{1, i}(1 : 2, :);
    B = DX_tmp{2, i}(1 : 2, :);

    c(i) = min([sum(sum((A - B) .^ 2) .^ 0.5),...
        sum(sum((diag([-1, 1]) * A - B) .^ 2) .^ 0.5),...
        sum(sum((diag([1, -1]) * A - B) .^ 2) .^ 0.5),...
        sum(sum((diag([-1, -1]) * A - B) .^ 2) .^ 0.5)]);
end

figure, semilogx(abs(c) / abs(c(2)), 'linewidth', 3);

xlabel('Log(#NNs)')
ylabel('L2 norm')
%yticks([0 1])