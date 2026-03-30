% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function a = angleVect(v1,v2,v3)
%
%       a = angleVect(v1,v2)
%
% return the angle (v1,v2)
% Inputs
%   - v1,v2 two vectors of size (3x1)
%   - v3 fixes the orientation to get the sign
%
% Outputs
%   - a the angle

v1 = normc(v1);
v2 = normc(v2);
v3 = normc(v3);

a = acos(v1'*v2)*sign(sin(acos(v1'*v2))*(cross(v1,v2)'*v3));

