% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function [d,al] = newPaperComputeD(pp,sh)
% compute the distance between the rules of the parameter set pp
% Inputs:
%   - pp: a structure giving the paper model parameter
%   - rg: a n-vector giving the region of the rules
%
% outputs:
%   - d: a n-vector giving the distance between the rules
%   - al: a (4xn) matrix giving the arc length of the previous and next points [p1;p2;n1;n2]

n = length(pp.be);

% add points between the rules
als = sort(pp.al(:)');
altot = pchip(-1:2*n+2,[als([end-1 end])-1 als als([1 2])+1],.5:.5:2*n+.5);

d = zeros(1,n);
al = zeros(4,n);
% compute the points and the distances
for i=1:n
    alpn = [altot(find(pp.al(1,i)>altot,1,'last')) altot(find(pp.al(1,i)<altot,1,'first')); ...
            altot(find(pp.al(2,i)<altot,1,'first')) altot(find(pp.al(2,i)>altot,1,'last'))];
    
    m = al2m(alpn([1 2 4 3]),sh);
    d(i) = polyarea(m(1,:),m(2,:))/mean([norm(m(:,2)-m(:,1)) norm(m(:,4)-m(:,3))]);
    al(:,i) = alpn(:);
end
