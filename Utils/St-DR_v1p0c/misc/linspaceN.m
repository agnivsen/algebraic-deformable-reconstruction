% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function X = linspaceN(X1 , X2 , varargin)
%
%   X = linspaceN(X1 , X2 , varargin)
%
% generates a row vector of 100 linearly
% equally spaced matrix between X1 and X2.


if nargin == 0;
    n = 100;
else
    n = varargin{1};
end
 n_elts = prod(size(X1));

 for i=1:n_elts
     X(i,:) = linspace(X1(i) , X2(i) , n);
 end
 

X = squeeze(reshape(X , [size(X1) n]));

if size(X,2) == 1 &&  ndims(X) == 2
    X = X';
end