% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function moveDelete

global figure_model rulingPositions rulingHandles ind_to_delete

set(figure_model,'pointer','crosshair')

th = .001;

figure(figure_model)
m = get(gca,'currentpoint');
m = m(1,1:2)';

% distance computation
n = size(rulingPositions,2);
d = [sum((rulingPositions(1:2,:) - repmat(m,1,n)).^2) ; ...
     sum((rulingPositions(3:4,:) - repmat(m,1,n)).^2) ; ...
     sum((rulingPositions(5:6,:) - repmat(m,1,n)).^2) ; ...
     th*ones(1,n)];
 
m = min(d(:));

set(rulingHandles,'color','k');
if m<th
    [i,ind_to_delete] = find(d==m);
    set(rulingHandles(:,ind_to_delete),'color','r');
end

