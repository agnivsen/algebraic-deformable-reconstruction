% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function [paper, rigidTransformation, points2D, corners2D] = paperRefinement(paperInit,sequenceIn,rigidTransformationInit, points2DInit, corners2DInit)

global sequence cameraMatrix
sequence = sequenceIn;

n_points = size(points2DInit,2);
ppInit = paperInit.ppnew;
n_parameters = length(ppInit.be);
n_borders = 25;
n_extra = 1;
cornerWeight = 1;

% form the parameter vector
% rigid motion
P(1:3) = rigidTransformationInit(:,4);
[u,d,v] = svd(rigidTransformationInit(:,1:3));
P(4) = d(1);
P(5:7) = find_angles(u*v');

% paper model
P(8) = ppInit.r;
% get the new parameterization to avoid crossing
ppInitValid = paperParameterisationConversion(ppInit,'valid_al_be');

P(9:9+2*n_parameters) = ppInitValid.ald(:);
P(10+2*n_parameters:9+3*n_parameters) = ppInitValid.be(:);

% points parameterisation
P(10+3*n_parameters:9+3*n_parameters+2*n_points) = points2DInit(:);

% jacobian matrix pattern
n_reprojectionPerPoint = sum(sequence.measurementMatrixPattern);
cs = cumsum(sequence.measurementMatrixPattern);
iind = cs + repmat([0 cumsum(n_reprojectionPerPoint(1:end-1))],size(sequence.measurementMatrixPattern,1),1);
iind = iind.*sequence.measurementMatrixPattern;
jind = repmat(1:n_points,size(sequence.measurementMatrixPattern,1),1).*sequence.measurementMatrixPattern;

J = zeros(sum(n_reprojectionPerPoint),2*n_points);
J(sub2ind(size(J),iind(iind~=0),2*jind(iind~=0)-1)) = 1;
J(sub2ind(size(J),iind(iind~=0),2*jind(iind~=0))) = 1;

n_reprojectionCorners = length(cat(2,sequence.corners.frameIndices))*2;

J = [ones(size(J,1),9+3*n_parameters) J ; ...
    ones(n_reprojectionCorners,9+3*n_parameters) zeros(n_reprojectionCorners,2*n_points)];

% non linear minimisation
sequence = sequenceIn;
cameraMatrix = cat(1,sequence.cameras.mat);

% compute the boujou rms error
n_cameras = size(cameraMatrix,1)/3;
ind_xy = vect([ones(2,n_cameras) ; zeros(1,n_cameras)])==1;
ind_z = ~ind_xy;
projectedPoints = (cameraMatrix(ind_xy,:) * [sequence.points3D ; ones(1,n_points)])./...
    (reshape(repmat(vect(cameraMatrix(ind_z,:))',2,1),2*n_cameras,4) * [sequence.points3D ; ones(1,n_points)]);

e = projectedPoints(sequence.measurementMatrixPattern==1)...
    - sequence.measurementMatrix(sequence.measurementMatrixPattern==1);
% add the error for the corners
corners3D = cat(2,sequence.corners.xyz);
for i=1:4
    n_cameras_i = length(sequence.corners(i).frameIndices);
    cameraMatrix_i = cat(1,sequence.cameras(sequence.corners(i).frameIndices).mat);
    ind_xy_i = vect([ones(2,n_cameras_i) ; zeros(1,n_cameras_i)])==1;
    ind_z_i = ~ind_xy_i;
    projectedCorner_i = (cameraMatrix_i(ind_xy_i,:) * [corners3D(:,i) ; 1])./...
        (reshape(repmat(vect(cameraMatrix_i(ind_z_i,:))',2,1),2*n_cameras_i,4) * [corners3D(:,i) ; 1]);
    e = [e ; projectedCorner_i - sequence.corners(i).uv(:)];
end
fprintf('\n*****\n')
fprintf('estimation with %d generatrices\n',n_parameters)
fprintf('*****\n')

fprintf('\n*****\n')
fprintf('boujou reconstruction rms error=%.2fpx\n',RMS(e));
fprintf('*****\n')

% evaluate the rms error before bundle adjustment
[nth,e] = errorPaperBS(P,ppInitValid.ind,n_parameters,n_points,n_extra,cornerWeight);
fprintf('\n*****\n')
fprintf('rms error before bundle adjustment=%.2fpx\n',RMS(e));
fprintf('*****\n')

% do the bundle adjustment
option = optimset('display' , 'iter', 'LargeScale', 'on', 'MaxIter', 50 , 'jacobPattern',J);
P = lsqnonlin(@(P)errorPaperBS(P,ppInitValid.ind,n_parameters,n_points,0,cornerWeight),P,[],[],option);
if n_extra~=0
    P = lsqnonlin(@(P)errorPaperBS(P,ppInitValid.ind,n_parameters,n_points,n_extra,cornerWeight),P,[],[],option);
end
[nth,e,paper, rigidTransformation, points2D, corners2D] = ...
    errorPaperBS(P,ppInitValid.ind,n_parameters,n_points,n_extra,cornerWeight);

fprintf('\n*****\n')
fprintf('rms error after bundle adjustment=%.2fpx\n',RMS(e));
fprintf('*****\n')

% plot the 3D figure
nx = round(n_borders*sqrt(paper.ppnew.r));
ny = round(n_borders*sqrt(1/paper.ppnew.r));
Minit = paperMesh(paperInit,nx,ny);
Minit = translate_mesh(rotate_mesh(Minit,rigidTransformationInit(:,1:3)),rigidTransformationInit(:,4));
M = paperMesh(paper,nx,ny);
M = translate_mesh(rotate_mesh(M,rigidTransformation(:,1:3)),rigidTransformation(:,4));

% compute the rulings
m1Init = rigidTransformationInit(:,1:3)*uv2M(paperInit.uv_rules(1:2,:),paperInit) + repmat(rigidTransformationInit(:,4),1,size(paperInit.uv_rules,2));
m2Init = rigidTransformationInit(:,1:3)*uv2M(paperInit.uv_rules(3:4,:),paperInit) + repmat(rigidTransformationInit(:,4),1,size(paperInit.uv_rules,2));
m1 = rigidTransformation(:,1:3)*uv2M(paper.uv_rules(1:2,:),paper) + repmat(rigidTransformation(:,4),1,size(paper.uv_rules,2));
m2 = rigidTransformation(:,1:3)*uv2M(paper.uv_rules(3:4,:),paper) + repmat(rigidTransformation(:,4),1,size(paper.uv_rules,2));

% uncomment to view the 3D meshes
% figure
% plotMesh(Minit,'k')
% hold on
% plotMesh(M,'g')
% plot3(sequence.points3D(1,:),sequence.points3D(2,:),sequence.points3D(3,:),'b.','markersize',10)
% plot3(corners3D(1,:),corners3D(2,:),corners3D(3,:),'r.','markersize',15)
% text(corners3D(1,:),corners3D(2,:),corners3D(3,:),num2str((1:4)'))
% plot3([m1Init(1,:) ; m2Init(1,:)],[m1Init(2,:) ; m2Init(2,:)],[m1Init(3,:) ; m2Init(3,:)],'c','linewidth',3)
% plot3([m1(1,:) ; m2(1,:)],[m1(2,:) ; m2(2,:)],[m1(3,:) ; m2(3,:)],'r','linewidth',3)
% end uncomment

for i=[1 length(sequence.frameList)]
    frame_i = imread([sequence.path 'frames/' sequence.frameList{i}]);
    figure
    imshow(frame_i)
    hold on
    plotMesh(projectMesh(M,sequence.cameras(i).mat),'w','linewidth',2)
end



