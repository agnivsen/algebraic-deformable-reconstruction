% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function [eWW,e,paper, rigidTransformation, points2D, corners2D] = errorPaperBS(P,indValid,n_parameters,n_points,n_extra,cornerWeight)

global cameraMatrix sequence

% extract the parameter set from the parameter vector
rigidTransformation = [P(4)*rotation(P(5:7)) P(1:3)'];


% paper model
ppValid.r = P(8);
ppValid.ald = P(9:9+2*n_parameters);
ppValid.be = P(10+2*n_parameters:9+3*n_parameters);
ppValid.ind = indValid;
ppValid.type = 'valid_al_be';

% points parameterisation
points2D = reshape(P(10+3*n_parameters:9+3*n_parameters+2*n_points),2,n_points);

% construct the paper
pp = paperParameterisationConversion(ppValid,'al_be',n_extra);
paper = newPaperFast(pp);

% compute the 3D points
points3D = rigidTransformation(:,1:3)*uv2M(points2D,paper) + repmat(rigidTransformation(:,4),1,n_points);

% project the 3D points
n_cameras = size(cameraMatrix,1)/3;
ind_xy = vect([ones(2,n_cameras) ; zeros(1,n_cameras)])==1;
ind_z = ~ind_xy;
projectedPoints = (cameraMatrix(ind_xy,:) * [points3D ; ones(1,n_points)])./...
    (reshape(repmat(vect(cameraMatrix(ind_z,:))',2,1),2*n_cameras,4) * [points3D ; ones(1,n_points)]);

% form the error vector
e = projectedPoints(sequence.measurementMatrixPattern==1)...
    - sequence.measurementMatrix(sequence.measurementMatrixPattern==1);

% project the corners
le = length(e);
corners2D = paper.sh.c;
corners3D = rigidTransformation(:,1:3)*uv2M(corners2D,paper) + repmat(rigidTransformation(:,4),1,4);
for i=1:4
    n_cameras_i = length(sequence.corners(i).frameIndices);
    cameraMatrix_i = cat(1,sequence.cameras(sequence.corners(i).frameIndices).mat);
    ind_xy_i = vect([ones(2,n_cameras_i) ; zeros(1,n_cameras_i)])==1;
    ind_z_i = ~ind_xy_i;
    projectedCorner_i = (cameraMatrix_i(ind_xy_i,:) * [corners3D(:,i) ; 1])./...
    (reshape(repmat(vect(cameraMatrix_i(ind_z_i,:))',2,1),2*n_cameras_i,4) * [corners3D(:,i) ; 1]);
    e = [e ; projectedCorner_i - sequence.corners(i).uv(:)];
end

eWW = e;
eWW(le+1:end) = cornerWeight*eWW(le+1:end);