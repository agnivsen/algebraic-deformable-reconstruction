% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function sg = signMP(M)
% return the signs of elements of M. The sign of 0 is 1.

sg = sign(M);
sg(sg==0) = 1;
