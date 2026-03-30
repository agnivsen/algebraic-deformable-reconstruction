close all; clear all; clear vars;

% Add path to CVX
addpath('alg-SfT/');
addpath('Data/');
addpath(genpath('Utils/'));

nNeighbors = 12; nFiles = 5; gridSize = [5 6]; noise = 0.0; showDebugPlots = false; equiarealityStrength = 0.0;
[Data, nng, K, ~, template] =  syntheticDataGenerator(nFiles, gridSize, nNeighbors, noise, showDebugPlots, equiarealityStrength);

deformationModel = 'Equiareal'; % Options are: 'Conformal' and 'Equiareal'

forwardDepthsWeights = 19; 
parameterSet = {10, 1, 0.1, forwardDepthsWeights}; % isoWt, traceWt, mdhWt (MOD), fdWt


huberThreshold = 0.1; precisionVal = 1.0e-03;


f = figure(1);
errList = [];

for index = 1:nFiles

    fprintf('<strong>HULK</strong> dataset: reconstructing image <strong>%d</strong> of %d\n\n', index, nFiles);

    eps = 0; solver = 'mosek'; mode = 'approximate';  

    if(strcmp(deformationModel, 'Conformal'))
        [reconstruction, solutionStatus] = algSfT_mConf(template, Data.p(index).p(1:2,:).', K, nng, eps, solver, parameterSet, mode, precisionVal);
    else
        [reconstruction, solutionStatus] = algSfT_mEqui(template, Data.p(index).p(1:2,:).', K, nng, eps, solver, parameterSet, mode, precisionVal);
    end

    [rms] = plotResults(reconstruction, Data.Pgth(index).P.');
    pause(0.1); hold off;

    errList = [errList rms];

end

fprintf('\n\nMean of RMSE (across all sequences) = <strong>%d</strong> mm\n\n', round(mean(errList)*1000, 2));
