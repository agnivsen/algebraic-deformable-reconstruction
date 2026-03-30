function [M,ind,J] = uv2M(uv,paper)
% [M,ind,J] = uv2M(uv,paper)
%
% Maps a set of points in uv coordinates to their 3D locations, and
% computes the embedding's Jacobian at these points.
%
% Inputs:
%
%  - uv [(2 x m) matrix]
%   The 2D points.
%
%  - paper [struct]
%   The paper model.
%
% Outputs:
%
%  - M [(3 x m) matrix]
%   The 3D points.
%  
%  - ind [m-vector]
%   The index of the cells where the points lie, from the paper model.
%
%  - J [(6 x m) matrix]
%   The Jacobian matrices, vectorized column-wise (ie the first 3
%   components are the embedding's u partial derivatives).

% number of points to embed
nuv = size(uv,2);

% number of rulings in the model
n = length(paper.ppnew.be);

% find the cells where the points lie in
% ind is (1 x nuv) and gives the indices of the cell where each point lies
test = signMP(paper.test.A*uv + repmat(paper.test.b,1,nuv));
sc = paper.test.mask * test;
[~,ind] = max(sc);

% embed the points to the reference plane with zero depth and in
% homogeneous coordinates
uvh = [uv;zeros(1,nuv);ones(1,nuv)];

% note: paper.RT is a matrix holding the rigid transformation for each cell
% (a (3 x 4) matrix [R T]). These are placed along the diagonal (paper.RT 
% is a block-diagonal matrix of 3x4 blocks).
% the number of cells is n+1 (one more cell than the number of rulings), so
% paper.RT is (3*(n+1) x 4*(n+1)).

% compute the columns in the matrix paper.RT for each uv (depends on the
% cell they lie in, computed as ind).
ind1 = sub2ind([4*(n+1),nuv],linspaceN(4*(ind-1)+1,4*ind,4)',repmat(1:nuv,4,1));
% compute the rows in the matrix paper.RT for each uv (depends on the
% cell they lie in, computed as ind).
ind2 = sub2ind([3*(n+1),nuv],linspaceN(3*(ind-1)+1,3*ind,3)',repmat(1:nuv,3,1));

% move the uvh points to the right rows (given by ind1)
M = zeros(4*(n+1),nuv);
M(ind1) = uvh;

% apply the transformations
M = paper.RT * sparse(M);

% the results are in ind2
M = full(M(ind2));

% form the index in paper.RT of the top-left element of the transformation
% for each point
ind3 = sub2ind([3*(n+1),4*(n+1)],3*(ind-1)+1,4*(ind-1)+1);

% create the indices for the whole blocks
ind4 = [ ind3 ; ind3+1 ; ind3+2 ; ind3+3*(n+1) ; ind3+3*(n+1)+1 ; ind3+3*(n+1)+2 ];

% retrieve the Jacobians as the first two columns of the rotations
J = zeros(6,nuv);
J(:) = full(paper.RT(ind4(:)));


