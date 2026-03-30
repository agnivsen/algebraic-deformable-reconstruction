% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function mesh2D = projectMesh(mesh3D , P)
% compute the projection of a mesh by the projective camera matrix P

Q = reshape(mesh3D , prod(size(mesh3D))/3 , 3);
q = projectPoints(Q' , P)';

mesh2D = reshape(q , [size(mesh3D , 1) , size(mesh3D , 2) , 2]);