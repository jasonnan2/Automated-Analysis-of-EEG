function fdr_mat = fdr_matCorrect(p_mat)
% Preforms fdr correction on upper triangular portion of a symmetric matrix

% Assume pmat is your symmetric p-value matrix (NxN)
N = size(p_mat, 1);

% Extract upper triangular indices (including diagonal)
idx = find(triu(ones(N), 0));

% Get the p-values in the upper triangle
pvals = p_mat(idx);

% Apply FDR correction (e.g., Benjamini-Hochberg)
adj_pvals= fdr(pvals);

% Now write corrected values back into a matrix
fdr_mat = zeros(N);           % initialize with NaNs if desired
fdr_mat(idx) = adj_pvals;   % fill upper triangle

% Optionally mirror it to make it symmetric
fdr_mat = fdr_mat + triu(fdr_mat,1)';% mirror upper to lower
end