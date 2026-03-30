% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function R=rotation(angle);
%
% R=rotation(angle);
% R is the rotation matrix defined by the 3 angles

    R=[cos(angle(1))*cos(angle(2)),cos(angle(1))*sin(angle(2))*sin(angle(3))-sin(angle(1))*cos(angle(3)),cos(angle(1))*sin(angle(2))*cos(angle(3))+sin(angle(1))*sin(angle(3));...
       sin(angle(1))*cos(angle(2)),sin(angle(1))*sin(angle(2))*sin(angle(3))+cos(angle(1))*cos(angle(3)),sin(angle(1))*sin(angle(2))*cos(angle(3))-cos(angle(1))*sin(angle(3));...
       -sin(angle(2)),cos(angle(2))*sin(angle(3)),cos(angle(2))*cos(angle(3))];
