% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function e = RMS(A,B)
% Compute the RMS error between A and B
if ~exist('B' , 'var') B = 0; end
e = sqrt(mean((A(:)-B(:)).^2));