% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function demo_model
% this function demonstrates the paper model.

close all

% set the path
pathOld = path;
path([pwd '/model'],path);
path([pwd '/misc'],path)
path([pwd '/other'],path)


% generate a paper parameter structure
pp.type = 'al_be'; % the id of the parameterization
pp.r = 1.4; % width length ratio
pp.al = [0.57 0.6374 0.9118 ; 1.1469 0.8737 1.111]; % arc lengths of the rulings
pp.be = [-5 -7 4]; % normalized bending angles

% interpolate the rulings
pp = paperParameterisationConversion(pp,'al_be',5);

% compute the paper
paper = newPaperFast(pp);

% display the flat paper
plotFlatPaper(paper);

% evaluate the 3D mesh
mesh3D = paperMesh(paper,20);
figure
plotMesh(mesh3D,'k')
axis equal

% restore the path
path(pathOld)