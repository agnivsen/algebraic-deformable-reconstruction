% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function demo_randModel
% this function demonstrates the paper model.

close all

% set the path
pathOld = path;
path([pwd '/model'],path);
path([pwd '/misc'],path)
path([pwd '/other'],path)

% compute the paper
paper = randPaper;

% display the flat paper
plotFlatPaper(paper);

% evaluate the 3D mesh
mesh3D = paperMesh(paper,20);
figure
plotMesh(mesh3D,'k')
axis equal

% restore the path
path(pathOld)