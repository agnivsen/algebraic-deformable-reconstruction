% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function [newM,T] = meshRigidTransformation(M,Mref,s)

% s: scale flag, is 1 a scale transformation is evaluated
if ~exist('s','var')
    s = 1;
end

% function [newM,T,scl,ang,trl] = meshRigidTransformation(M,Mref)
sz = size(M);
n = sz(1)*sz(2);
Q = reshape(M,n,3)';
Qref = reshape(Mref,n,3)';
centerQ = mean(Q,2);
centerQref = mean(Qref,2);

Q_c = Q - repmat(centerQ,1,n);
Qref_c = Qref - repmat(centerQref,1,n);

if s
    v = mean(sqrt(sum(Q_c.^2)));
    vref = mean(sqrt(sum(Qref_c.^2)));

    
else
    vref = 1;
    v = 1;
end

Q_cs = vref/v * Q_c;

q = RotationUsingQuat(Q_cs, Qref_c, n);
R = quat2rot(q);
T = [eye(3) centerQref]*[R*vref/v*[eye(3) -centerQ] ; [0 0 0 1]];

newM = reshape((T*[Q;ones(1,n)])',[sz(1),sz(2),3]);

