% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function R = quat2rot(q);

% s = q(1);
% l = q(2);
% m = q(3);
% n = q(4);
% 
% R = [
%     s*s+l*l-m*m-n*n 2*(l*m-s*n) 2*(l*n+s*m)
%     2*(l*m+s*n) s*s-l*l+m*m-n*n 2*(m*n-s*l)
%     2*(l*n-s*m) 2*(m*n+s*l) s*s-l*l-m*m+n*n
% ];
% %return;

R = [
    1-2*(q(3)^2+q(4)^2), 2*(q(2)*q(3)-q(1)*q(4)), 2*(q(2)*q(4)+q(1)*q(3))
    2*(q(2)*q(3)+q(1)*q(4)), 1-2*(q(4)^2+q(2)^2), 2*(q(3)*q(4)-q(1)*q(2))
    2*(q(2)*q(4)-q(1)*q(3)), 2*(q(3)*q(4)+q(1)*q(2)), 1-2*(q(2)^2+q(3)^2)
];
