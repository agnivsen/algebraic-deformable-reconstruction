% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function ppNew = paperOptimalParameterisationWB(paper,th)

% increase the number of parameters to get a better approximation of the shape
if ~exist('th','var')
    th = .1;
end

% compute a new set of parameter based on normalized curve fitting
pp = paper.ppnew;
d = newPaperComputeD(pp,make_shape(pp.r));

alNew = [];
beNew = [];
for i=1:max(paper.rg)
    ind = paper.rg==i;
    [CP_i,altmp,betmp] = curveFitting(pp.al(1,ind),pp.al(2,ind),pp.be(ind).*d(ind),th);
    CP_i = unique(CP_i','rows')';
    alNew = [alNew,CP_i(1:2,:)];
    beNew = [beNew,CP_i(3,:)];
end

ppNew = pp;
ppNew.al = alNew;
ppNew.be = beNew;
[ppNew.al , ppNew.be] = uncross(ppNew.al , ppNew.be);
