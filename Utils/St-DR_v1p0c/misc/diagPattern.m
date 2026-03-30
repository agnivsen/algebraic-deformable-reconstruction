% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function M = diagPattern(n,nr,nc)
%
%   M = diagPattern(n,nr,nc)
%
% build a (nrxnc) block diagonal matrix with n blocks

M = eye(n);
M = repmat(M(:)',nr,1);
M = reshape(M,n*nr,n);
M = repmat(M,nc,1);
M = reshape(M,n*nr,n*nc);

