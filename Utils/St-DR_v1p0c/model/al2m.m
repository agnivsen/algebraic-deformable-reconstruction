% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function m = al2m(al , sh)
%
%   m = al2m(al , sh)
%
% convert arc lengths to 2D according to the shape.
% Inputs:
% - al: arc lengths of the points
% - sh: shape
% 
% Outputs:
% - m: the 2d points

al = mod(al , 1);
m = [ppval(sh.ppx , al); ppval(sh.ppy , al)];
