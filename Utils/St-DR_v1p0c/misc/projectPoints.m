% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function q = projectPoints(Q , P)

q = P * [Q ; ones(1 , size(Q , 2))];
q = [q(1,:) ./ q(3,:) ; q(2,:) ./ q(3,:)];