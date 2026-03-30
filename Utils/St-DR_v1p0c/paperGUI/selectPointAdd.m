% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function selectPointAdd

global current_ruling current_point figure_model fsPoint current_plot current_rulingHandle sh rulingPositions rulings rulingHandles

% delete the temporary plot
if ~isnull(ishandle(current_plot))
    delete(current_plot(ishandle(current_plot)));
end

switch fsPoint
    case 'f'
        current_ruling = current_point;        
        fsPoint = 's';
        current_rulingHandle = plot(current_point(1),current_point(2),'k.','markersize',25);
    case 's'
        set(figure_model,'WindowButtonMotionFcn','moveParameter')
        set(figure_model,'WindowButtonDownFcn','downParameter')
        set(figure_model,'WindowButtonUpFcn','upParameter')
        set(figure_model,'pointer','fleur')

        
        current_ruling = [current_ruling current_point];
        current_al = m2al(current_ruling , sh);
        [current_al,ind] = sort(current_al);
        current_ruling = current_ruling(:,ind);
        
        current_plot = [];
        current_rulingHandle(2) = plot(current_point(1),current_point(2),'k.','markersize',25);
        current_rulingHandle(1:2) = current_rulingHandle(ind);
        
        current_rulingHandle(3) = plot(mean(current_ruling(1,:)),mean(current_ruling(2,:)),'k.','markersize',25);
        
        

        
        rulingPositions = [rulingPositions [current_ruling(:);mean(current_ruling,2)]];        
        current_rulingHandle(4) = plot(current_ruling(1,ind),current_ruling(2,ind),'k');
        
        rulings = [rulings [current_al' ; 0]]; % ruling initialisation (bending angle is set to 0)
        rulingHandles = [rulingHandles current_rulingHandle'];
        
        plotPaper
end


