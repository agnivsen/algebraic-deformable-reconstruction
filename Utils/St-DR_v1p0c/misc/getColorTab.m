% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function co = getColorTab(n)
% co = getColorTab(n)
%
% return a tab nx3 where each row is a different rgb color

h = figure();
c = colormap; % get the default colormap
close(h)

ind = round(linspace(1,size(c,1),n));
co = c(ind,:);

if n> size(c,1)
    warning('There is not enough color to represent each value');
end