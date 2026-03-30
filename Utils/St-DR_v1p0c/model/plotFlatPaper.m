% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function h = plotFlatPaper(paper,delta,c)


if ~exist('c','var')
    c = 1;
end

if ~exist('delta','var')
    delta = 1;
end

if c==1
    if max(paper.rg)==3
        colors = [0 0 1;0 1 1;1 1 0];
    else
        colors = getColorTab(max(paper.rg));
    end
elseif c==-1
    colors = getColorTab(max(paper.rg));
    colors = colors([5 2 4 3 1],:);
end


h = figure;
hold on
for i=1:delta:size(paper.uv_rules,2)
    if c
        plot(paper.uv_rules([1 3],i),paper.uv_rules([2 4],i),'color',colors(paper.rg(i),:),'linewidth',3)
    else
        plot(paper.uv_rules([1 3],i),paper.uv_rules([2 4],i),'g')
    end
end

plot(paper.sh.poly(1,:),paper.sh.poly(2,:),'k','linewidth',3);
