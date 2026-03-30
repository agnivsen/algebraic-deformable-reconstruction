% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function [I,l] = getMultiInter(A,B,C,D)
% compute the point I as the intersection of AB and CD, such that I = A+l*(B-A)
% return nan if AB and CD are colinear
% A,B,C,D can be (2xn) matrices of points


v1 = A-D;
v2 = C-D;
v3 = B-A;

l1 = (v1(1,:).*v2(2,:)-v1(2,:).*v2(1,:));
l2 = (v2(1,:).*v3(2,:)-v2(2,:).*v3(1,:));

ind = l2==0;
l2(ind) = nan;
l = l1./l2;
I = A + repmat(l,2,1).*v3;