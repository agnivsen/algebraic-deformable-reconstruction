% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function ind_al = get_intersections(al)
% al is a 2xn matrix. each column represent a line, parameterized by its 2 arc lenghts.
% Inputs
% - al: the matrix of line coordinates
% Outputs:
% - ind: a 2xm matrix. m is the number of intersections. ind(:,k) are the indices line where the kth intersection occurs  
al = sort(mod(al,1));
[al(1,:),indo] = sort(al(1,:));
al(2,:) = al(2,indo);
%[nth,indo] = sort(indo);

th_base = eps;
th_intersect = th_base; % the threshold to detect the crossing rules
th_uncross = 2*th_base; % the threshold to avoid the crossing of rules. It must be >= to th_intersect
if th_intersect > th_uncross
    error('The thresholds to detect crossing and to correct crossing are not compatible');
end
th_mesh = 3*th_base; % the threshold to generate a mesh inside the paper

n = size(al , 2);

AL_ind = [nonzeros(tril(repmat([2:n]',1,n-1)))' ; nonzeros(tril(repmat([1:n-1],n-1,1)))'];
AL_comp = al(:,AL_ind(1,:));
AL_ref = repmat(al(2,AL_ind(2,:)),2,1);

AL_test = AL_comp - AL_ref;
ind = (AL_test(1,:)<th_intersect & AL_test(2,:)>-th_intersect);

ind_al = indo(AL_ind(:,ind));

% tricks to correct a bug
if size(ind_al,1)==1
    ind_al = ind_al';
end
