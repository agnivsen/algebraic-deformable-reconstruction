% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function R = rot_vect(v,phi)
% R = rot_vect(v,phi)
%
% Compute the rotation arround the vector v of the angle phi.
% Inputs:
% - v: a (3x1) matrix giving the axe vector of the rotation
% - phi: the angle of the rotation
%
% Output:
% - R: the (3x3) rotation matrix

v = normc(v(:));
q = [cos(phi/2), sin(phi/2)*v'];

R = [
    1-2*(q(3)^2+q(4)^2), 2*(q(2)*q(3)-q(1)*q(4)), 2*(q(2)*q(4)+q(1)*q(3))
    2*(q(2)*q(3)+q(1)*q(4)), 1-2*(q(4)^2+q(2)^2), 2*(q(3)*q(4)-q(1)*q(2))
    2*(q(2)*q(4)-q(1)*q(3)), 2*(q(3)*q(4)+q(1)*q(2)), 1-2*(q(2)^2+q(3)^2)
];