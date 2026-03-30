% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function movePointerAdd
global figure_model sh current_plot current_point fsPoint current_ruling



if ~isnull(ishandle(current_plot))
    delete(current_plot(ishandle(current_plot)));
    current_plot = [];
end
figure(figure_model)
m = get(gca,'currentpoint');
m = m(1,1:2)';

I = getMultiInter(sh.poly(:,1:4),sh.poly(:,2:5),zeros(2,4),repmat(m,1,4));
v = I'*m;
v(v<0) = inf;
[nth,ind] = min(v);
current_point = I(:,ind);
current_plot = [current_plot plot(current_point(1),current_point(2),'r.','markersize',25)];
if isequal(fsPoint,'s')
    current_plot = [current_plot plot([current_point(1) current_ruling(1)],[current_point(2) current_ruling(2)],'k--')];
end