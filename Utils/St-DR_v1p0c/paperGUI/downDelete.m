% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function downDelete

global figure_model ind_to_delete rulingHandles rulings rulingPositions


set(figure_model,'pointer','fleur')

delete(rulingHandles(:,ind_to_delete))
rulings(:,ind_to_delete) = [];
rulingHandles(:,ind_to_delete) = [];
rulingPositions(:,ind_to_delete) = [];

set(figure_model,'WindowButtonMotionFcn','moveParameter')
set(figure_model,'WindowButtonDownFcn','downParameter')
set(figure_model,'WindowButtonUpFcn','upParameter')
