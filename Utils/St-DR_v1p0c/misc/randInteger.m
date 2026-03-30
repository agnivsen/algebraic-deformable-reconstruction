% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function r = randInteger(a , b , varargin)

%
% r = randInteger(a , b , varargin)
% Uniformly distributed random integers in [a , b]
% See also rand
a = a - .5;
b = b + .5;
r = round((b - a) * rand(varargin{:}) + a * ones(varargin{:}));