% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function plotFlatPaperWithMeshInfo(paper,nv,f1,f2,cl)

tic
[M,m,ind] = paperMesh(paper,nv);
toc
size(m)
if exist('f1','var')
    figure(f1)
else
    figure
end

plotMesh(m,'y');
hold on
c = getColorTab(max(ind(:)));
for i=1:size(m,1)
    for j=1:size(m,2)
        text(m(i,j,1),m(i,j,2),num2str(ind(i,j)),'color',c(ind(i,j),:),'backgroundcolor',1-c(ind(i,j),:));
    end
end

plot(paper.sh.poly(1,:),paper.sh.poly(2,:),'b');
mo = al2m(0,paper.sh);
plot(mo(1),mo(2),'r.','markersize',10)

for i=1:size(paper.uv_rules,2)
    plot(paper.uv_rules([1 3],i),paper.uv_rules([2 4],i),'r')
    text(mean(paper.uv_rules([1 3],i)),mean(paper.uv_rules([2 4],i)),num2str(i),'color',[0 0 1],'backgroundcolor',[0 1 0])
end

if exist('f1','var')
    figure(f2)
else
    figure
end
if exist('cl','var')
    plotMesh(M,'color',cl)
else
    plotMesh(M,'k')
end

axis equal


