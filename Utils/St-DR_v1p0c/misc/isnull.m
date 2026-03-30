% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function n = isnull(A)
% n = isnull(A)
% 
% return 1 if all elements of A are null, 0 otherwise

n = isequal(A , zeros(size(A)));