% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

% threshold definitions for the geometric part of the paper model
th_base = eps;
th_intersect = th_base; % the threshold to detect the crossing rules
th_uncross = 2*th_base; % the threshold to avoid the crossing of rules. It must be >= to th_intersect
if th_intersect > th_uncross
    error('The thresholds to detect crossing and to correct crossing are not compatible');
end
th_mesh = 3*th_base; % the threshold to generate a mesh inside the paper