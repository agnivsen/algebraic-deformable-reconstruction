% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function mesh3D = translate_mesh(mesh3D , T)

vm = reshape(mesh3D , prod(size(mesh3D))/3 , 3)' + repmat(T , 1 , prod(size(mesh3D))/3);
mesh3D = reshape(vm' , size(mesh3D));