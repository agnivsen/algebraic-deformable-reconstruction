function [error, cloud1, cloud2] = RMSE3D(cloud1, cloud2)
%% RMSE3D : computes RMSE3D between two pointclouds
%
% PARAMETERS:
%   cloud1 :  first input pointcloud of dimension [N x 2] or [N x 3],
%                   although [N x M] dimensional pointcloud can be handled as well
%   cloud2 :  second input pointcloud of dimension [N x 2] or [N x 3],
%                   although [N x M] dimensional pointcloud can be handled as well
%
%   (The dimension of cloud1 and cloud2 needs to be same)
%
%-------------------------------------------------------------------------------------
%
% METHOD:
%   Does the following: sqrt{ (\sum_{i=1}^N ( (X_i1 - X_i2)^2 + (Y_i1 - Y_i2)^2 + (Z_i1 - Z_i2)^2 ) ) / N }
%
% RETURNS:
%   error: the computed RMSE


    assert(size(cloud1,1) == size(cloud2,1),'RMSE: the two pointcloud should be of same dimension');
    assert(size(cloud1,2) == size(cloud2,2),'RMSE: the two pointcloud should be of same dimension');
    
    delIndex = [];
    
    for iP = 1:size(cloud1,1)
        if( isnan(cloud1(iP,1)) || isnan(cloud1(iP,2)) || isnan(cloud1(iP,3)) || isinf(cloud1(iP,1)) || isinf(cloud1(iP,2)) || isinf(cloud1(iP,3)) || ...
                isnan(cloud2(iP,1)) || isnan(cloud2(iP,2)) || isnan(cloud2(iP,3)) || isinf(cloud2(iP,1)) || isinf(cloud2(iP,2)) || isinf(cloud2(iP,3))...
                || (abs(cloud1(iP,3)) < 10^(-20)) || (abs(cloud2(iP,3)) < 10^(-20)) )
            delIndex = [delIndex; iP];
        end
    end
    
    cloud1(delIndex,:) = [];
    cloud2(delIndex,:) = [];
    
    if(size(cloud1,2)>3)
        warning('RMSE: was expecting a point cloud of dimension [N x 2] or [N x 3]');
        fprintf('Instead, received pointcloud dimension = %d\n', size(cloud1,2));
    end

    squaredDiff = (cloud1.' - cloud2.').^2;
    sumOfSquaredDiff = (sum(squaredDiff));
    meanSumOfSquaredDiff = mean(sumOfSquaredDiff);
    rootMeanSumOfSquaredDiff = sqrt(meanSumOfSquaredDiff);
    error = rootMeanSumOfSquaredDiff;
end