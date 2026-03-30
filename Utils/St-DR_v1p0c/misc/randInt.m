% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function r = randInt(a , b , varargin)
%
% r = randInt(a , b , varargin)
% Uniformly distributed random numbers in [a , b]
% See also rand

r = (b - a) * rand(varargin{:}) + a * ones(varargin{:});
