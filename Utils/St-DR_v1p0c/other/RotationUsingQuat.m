% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function q = RotationUsingQuat(Q1, Q2, m);

%% Compute a rotation between the sets of m points Q1 and Q2, represented
%% by 3 x m matrices
%
% The sets of points must be centred on their centroid (translation
% compensated). The procedure thus minimizes the squared distance error.

A = zeros(4);
for j = 1:m
    X1 = Q1(:, j);
    X2 = Q2(:, j);
%     Aj = [
%       0 (X1-X2)'
%       X2-X1 mcrp(X2)+mcrp(X1)
%     ];
    Aj = [
      0 (X1-X2)'
      X2-X1 crMat(X2)+crMat(X1)
    ];
    A = A + Aj'*Aj;
end;
[u,d,v] = svd(A);
q = v(:, end);
