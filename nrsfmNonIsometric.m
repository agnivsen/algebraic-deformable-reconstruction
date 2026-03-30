close all; clear all; clear vars;

% Add path to CVX
addpath('alg-NRSfM/');
addpath('Data/');
addpath(genpath('Utils/'));

nNeighbors = 7; nFiles = 8; gridSize = [5 4]; noise = 0.0; showDebugPlots = false; equiarealityStrength = 0.0;
[Data, nng, K, ~, template] =  syntheticDataGenerator(nFiles, gridSize, nNeighbors, noise, showDebugPlots, equiarealityStrength);

deformationModel = 'Conformal'; % Options are: 'Conformal' and 'Equiareal'

forwardDepthsWeights = 19;
parameterSet = {10, 1, 0.1, forwardDepthsWeights}; % isoWt, traceWt, mdhWt (MOD), fdWt

precisionVal = 1.0e-03; lambda = 7; eps = 0.0; solver = 'mosek';


[reconsFD, ~] = NRSfM_IsometricZeroth(Data, nng, K, 'DeformationModel','mdh-extensible-isometric', ...
    'RescaleData', false, 'Solver', 'mosek');

if(strcmp(deformationModel, 'Conformal'))
    reconstruction = algNRSfM_mConf(Data, K, nng, reconsFD.reconstruction, reconsFD.summedGeodSquared, eps, solver, lambda);
else
    reconstruction = algNRSfM_mEqui(Data, K, nng, reconsFD.reconstruction, reconsFD.summedGeodSquared, eps, solver, lambda);
end

[errList] = plotResultsNRSfM(reconstruction, Data);

fprintf('\n\nMean of RMSE (across all sequences) = <strong>%d</strong> units\n\n', round(mean(errList)*1000, 2));