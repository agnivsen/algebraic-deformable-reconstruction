% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function [al , ba] = uncross(al , ba)
% change the arc length such that the rules defined by al don't cross each other
threshold_developable; % load the value of thresholds

ind = get_intersections(al);

while ~isempty(ind)
    ind = ind(:,1);
    M_al_1 = al(:,ind(1,:));
    M_al_2 = al(:,ind(2,:));
    M_al_1 = repmat(M_al_1 , 2 , 1);
    M_al_2 = [M_al_2(1,:) ; M_al_2(1,:) ; M_al_2(2,:) ; M_al_2(2,:)];
    d_al = abs(M_al_1 - M_al_2);

    [i,k] = find(d_al == min(min(d_al)));
    [i,j] = ind2sub([2,2],i);
    al1 = al(i,ind(1,k));
    al2 = al(j,ind(2,k));
    s = sign(al2-al1);
    s(s==0)=1;

    al(i,ind(1,k)) = mean([al1 al2])+s*th_uncross/2;
    al(j,ind(2,k)) = mean([al1 al2])-s*th_uncross/2;

    al = sort(al);
    [nth,ind] = sort(al(1,:));
    al = al(:,ind);
    ba = ba(ind);
    ind = get_intersections(al);

end