% Computes the effective rank of a matrix, ignoring
% singular values below a 'threshold' (default: 10^-11), which
% usually crops up due to numerical issues
function  [rank] = effectiveRank(A, threshold)
    if ~exist('threshold', 'var')
        threshold = 10^(-11);
    end
    s = svd(A);
    rank = sum(s/s(1) > threshold);
end