% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function closeAll_guiPaper

global figure_3D figure_model
if ishandle(figure_3D)
    close(figure_3D)
end
if ishandle(figure_model)
    close(figure_model)
end

% restore the path
load pathTmp
path(pathOld)