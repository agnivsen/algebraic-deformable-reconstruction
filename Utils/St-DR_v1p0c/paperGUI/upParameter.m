% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function upParameter

global isDown
global figure_model rulingPositions rulingHandles rulings current_plot current_bending coef_bending ind_changing current_ruling sh

isDown = 0;

if ~isnull(ishandle(current_plot))
    delete(current_plot(ishandle(current_plot)));
    current_plot = [];
end

if ~isempty(ind_changing)
    if ind_changing(1)==3
        rulingPositions(5:6,ind_changing(2)) = current_bending;
        delete(rulingHandles(ind_changing(1),ind_changing(2))); 
        rulingHandles(ind_changing(1),ind_changing(2)) = plot(current_bending(1),current_bending(2),'g.','markersize',25);
        rulings(3,ind_changing(2)) = 2*pi*(coef_bending-.5);
        
        ind_changing = [];
        coef_bending = [];
        current_bending = [];
    else
        % order the points
        current_al = m2al(current_ruling , sh);
        [current_al,ind] = sort(current_al);
        current_ruling = current_ruling(:,ind);

        % delete and plot the ruling
        current_plot = [];
        delete(rulingHandles(:,ind_changing(2)));
        current_rulingHandle(1) = plot(current_ruling(1,1),current_ruling(2,1),'k.','markersize',25);
        current_rulingHandle(2) = plot(current_ruling(1,2),current_ruling(2,2),'k.','markersize',25);
        current_bending = (rulings(3,ind_changing(2))/(2*pi) + .5)*diff(current_ruling,1,2) + current_ruling(:,1);
        current_rulingHandle(3) = plot(current_bending(1),current_bending(2),'k.','markersize',25);
        current_rulingHandle(4) = plot(current_ruling(1,ind),current_ruling(2,ind),'k');

        % save the new parameters
        rulingHandles(:,ind_changing(2)) = current_rulingHandle;
        rulingPositions(:,ind_changing(2)) = [current_ruling(:);current_bending];
        
        rulings(1:2,ind_changing(2)) = current_al';
        
        ind_changing = [];
        
    end
end
plotPaper

return
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
    set(rulingHandles(d==m),'color','g');
    set(figure_model,'pointer','circle')
else
    set(figure_model,'pointer','fleur')
end


