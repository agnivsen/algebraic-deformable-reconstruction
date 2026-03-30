% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function demo_modelGUI
% this function demonstrates the paper model by launching the GUI.

% set the path

pathOld = path;
save pathTmp pathOld
path([pwd '/model'],path);
path([pwd '/paperGUI'],path)
path([pwd '/misc'],path)
path([pwd '/other'],path)

% start the Graphical User Interface
paper_GUI

