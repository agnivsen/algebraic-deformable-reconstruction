% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function [ M ] = crMat(v);
% [ M ] = crMat(v);
%
% Returns the (3x3) matrix M such that Mz = v x z for all z \in \R^3.

M = [
    0 -v(3) v(2)
    v(3) 0 -v(1)
    -v(2) v(1) 0
];
