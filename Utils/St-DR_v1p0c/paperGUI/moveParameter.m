% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function moveParameter

global figure_model rulingPositions rulingHandles isDown current_bending current_plot ind_changing coef_bending sh current_ruling

if ~isempty(rulingPositions)
if ~isnull(ishandle(current_plot))
    delete(current_plot(ishandle(current_plot)));
    current_plot = [];
end


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
if md<th||isDown
    set(figure_model,'pointer','circle')    
    if isDown && ~isempty(ind_changing)
        if ind_changing(1)==3 % change the bending angle
            
            set(rulingHandles(ind_changing(1),ind_changing(2)),'color',[.6 .6 .6]);
            v = rulingPositions(3:4,ind_changing(2))-rulingPositions(1:2,ind_changing(2));
            coef_bending = (m-rulingPositions(1:2,ind_changing(2)))'*v/ norm(v)^2;
            if coef_bending<0
                coef_bending = 0;
            elseif coef_bending>1
                coef_bending = 1;
            end
            current_bending = rulingPositions(1:2,ind_changing(2)) + coef_bending * v;
            current_plot = [current_plot plot(current_bending(1),current_bending(2),'r.','markersize',25)];
            
            % plot ticks along the rulings to give a scale idea about the bending angle
            current_tick = repmat(rulingPositions(1:2,ind_changing(2)),1,3) + [.25 .5 .75 ; .25 .5 .75] .* repmat(v,1,3);
            current_plot = [current_plot text(current_tick(1,:),current_tick(2,:),{'-' '0' '+'})'];
            
        else % change the ruling position
            set(rulingHandles(:,ind_changing(2)),'color',[.6 .6 .6]);
            I = getMultiInter(sh.poly(:,1:4),sh.poly(:,2:5),zeros(2,4),repmat(m,1,4));
            v = I'*m;
            v(v<0) = inf;
            [nth,ind] = min(v);
            current_point = I(:,ind);
            current_ruling = [current_point rulingPositions(2*(3-ind_changing(1))-1:2*(3-ind_changing(1)),ind_changing(2))];
            current_plot = [current_plot plot(current_point(1),current_point(2),'r.','markersize',25)];
            current_plot = [current_plot plot(current_ruling(1,:),current_ruling(2,:),'k--')];
        end
    else
        set(rulingHandles(d==md),'color','g');
    end    
else
    isDown = 0;
    set(figure_model,'pointer','fleur')
end


end