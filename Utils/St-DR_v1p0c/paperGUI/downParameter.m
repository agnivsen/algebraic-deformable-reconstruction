% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function downParameter


global figure_model rulingPositions rulingHandles ind_changing
global isDown
isDown = 1;


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
 
md = min(d(:));


set(rulingHandles,'color','k');
if md<th
    [ind_changing(1),ind_changing(2)] = find(d==md);    
    set(rulingHandles(d==md),'color','r');
else
    set(figure_model,'pointer','fleur')
end


