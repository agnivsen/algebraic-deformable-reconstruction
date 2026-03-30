% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function demo_estimatePaper
% this function demonstrates the estimation algorithm
close all

addpath('/srv/public/asengupta/Code/SfT_Baseline/GTPSW_v1p2/GTPSW_v1p2/');

% set the path
pathOld = path;
path([pwd '/model'],path);
path([pwd '/estimation'],path)
path([pwd '/misc'],path)
path([pwd '/other'],path)

% load the sequence (the planar parameteization has been previously computed
load example/data.mat
sequence.path = [pwd '/example/'];

% initialize the model
[paperInit, rigidTransformation, points2D, corners2D] = paperInitialization(sequence.points3D,[sequence.corners.xyz],planarParameterization.points2D,planarParameterization.corners2D);

% refine the model
[paperInit2, rigidTransformation2, points2D2, corners2D2] = paperRefinement(paperInit,sequence, rigidTransformation, points2D, corners2D);

% restore the path
path(pathOld)