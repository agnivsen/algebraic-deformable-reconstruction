% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function [pp,rg] = oabtr(pp)
%   [pp,rg] = oabtr(pp)
% order and bring together regions
% for a al/be parameterisation

% order the arc lengths
[pp.al,ind] = sort(mod(pp.al,1));
sgn = 2*(1.5-ind(1,:));
% change the sign of bending if needed
pp.be = pp.be.*sgn;
% order the rules
[pp.al(1,:),ind] = sort(pp.al(1,:));
pp.al(2,:) = pp.al(2,ind);
pp.be = pp.be(ind);

% reorder the rules
[nth,ind] = sort(pp.al(:));
ind_first = ind(find(diff(ind) == 1,1,'first')+1)/2;
pp.al(1,1:ind_first) = pp.al(1,1:ind_first)+1;

% order the arc lengths
[pp.al,ind] = sort(pp.al);
sgn = 2*(1.5-ind(1,:));
% change the sign of bending if needed
pp.be = pp.be.*sgn;

% order the rules
[pp.al(1,:),ind] = sort(pp.al(1,:));
pp.al(2,:) = pp.al(2,ind);
pp.be = pp.be(ind);

if nargout>1
    % compute regions
    [nth, ind] = sort(pp.al(2,:),'descend');
    rg = cumsum([1 diff(ind)~=1]);
    rg(ind) = rg;
    rg = cumsum([1 diff(rg)~=0]);
end

pp.type = 'al_be';