% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function [paper, rigidTransformation, points2D, corners2D, interp] = paperInitialization(points3D,corners3D,points2D,corners2D)
% initialize the paper model from a 3D poins cloud.
% Input:
%  -points3D: a (3xn) matrix giving the 3D points
%  -corners3D: a (3x4) matrix giving the 3D corners
%  -points2D: a (2xn) matrix giving the 3D points
%  -corners2D: a (2x4) matrix giving the 3D corners
%  
% Output:
%  -paper: the paper model
%  -rigidTransformation: the rigid transformation to apply to the model to fit it to the 3D points
%  -points2D: a (2xnout) matrix giving the planar parameterization of the points on the paper (nout can be not 
% equal to n since points that are flaged out of the shape are removed).
%  -corners2D: a (2x4) matrix giving the 3D corners
%
% See also newPaperFast,


%%
%% Settings

n_borders = 25; % the number of points along each side of the boundary to compute the rulings
n_samples = 50; % the number of samples along a ruling to evaluate its score
threshold_closeness = .05; % the threshold that define the closeness between rulings
n_interpolation = 128; % the size of the ruling position graph
n_rulingCurveMin = 10; % minimal number of rulings to form a curve
n_maxRulings = 10; % the maximal number of rulings starting from one point in the detection (filtering)
th_curveFitting = .02; % maximal error allowed in the curve fitting. Smaller value means more control points.
TPSregularization = 5e-4; % the TPS regularization


%%
%% 3D shape computation

%% Ruling detection
points2D_all = [points2D corners2D];
points3D_all = [points3D corners3D];

TPS = points3D_all*GTPSW_TPS_compute_FDPMat(points2D_all',[],TPSregularization)';

% get the points along the border
dx = 2*corners2D(1,2);
dy = 2*corners2D(2,3);

nx = round(n_borders*sqrt(dx/dy));
ny = round(n_borders/sqrt(dx/dy));
px = linspace(corners2D(1,1),corners2D(1,2),nx);
py = linspace(corners2D(2,2),corners2D(2,3),ny);
borderPoints = [[px(1:end-1);corners2D(2,2)*ones(1,nx-1)] ...
                 [corners2D(1,2)*ones(1,ny-1);py(1:end-1)] ...
                 [px(end:-1:2);corners2D(2,3)*ones(1,nx-1)] ...
                 [corners2D(1,1)*ones(1,ny-1);py(end:-1:2)]];

n_bordersPoints = size(borderPoints,2);
rulingScore = [1:n_bordersPoints ; nan(1,n_bordersPoints) ; inf(1,n_bordersPoints)]; % the two first 
% lines are ruling point indices and the third one is the score

% build the correspondence matrix (border point indice correspondances for candidate rulings)
candidateRulings = [[ones(1,nx+ny-1) ; nx:2*nx+ny-2]...
                    [vect(repmat(2:nx-1,nx+2*ny-4,1))' ; repmat(nx+1:2*nx+2*ny-4,1,nx-2)]...
                    [repmat(nx,1,nx+ny-2) ; nx+ny-1:2*nx+2*ny-4]...
                    [vect(repmat(nx+1:nx+ny-2,nx+ny-3,1))' ; repmat(nx+ny:2*nx+2*ny-4,1,ny-2)]...
                    [repmat(nx+ny-1,1,ny-1) ; 2*nx+ny-2:2*nx+2*ny-4]...
                    [vect(repmat(nx+ny:2*nx+ny-3,ny-2,1))' ; repmat(2*nx+ny-1:2*nx+2*ny-4,1,nx-2)]];

% compute the score of the candidate ruling and update the rulingScore matrix
ls = [linspace(1,0,n_samples);linspace(0,1,n_samples)]; 

% to show a progress bar
hwb = waitbar(0,'Ruling Detection Progress');
hAxes = findobj(hwb,'type','axes');
hTitle = get(hAxes,'title');
n_candidateRulings = size(candidateRulings,2);

for i=1:n_candidateRulings
    % compute the rulingVectors and the surface normal vectors
    ruling2D = borderPoints(:,candidateRulings(:,i))*ls;
%    [W,Wu,Wv] = GTPSW_TPS_kernelize_fafd(points2D_all',ruling2D');
    W = GTPSW_TPS_kernelize(points2D_all',ruling2D');
    [Wu,Wv] = GTPSW_TPS_kernelize_der1(points2D_all',ruling2D');
    
    ruling3D = TPS*W';
    rulingVectors = normc(diff(ruling3D,1,2));
    normalVectors = normc(cross(TPS*Wu',TPS*Wv'));

    % evaluate the scores
    scA = var(acos(normc(mean(rulingVectors,2))'*rulingVectors)); % the points are aligned
    scP = var(acos(normc(mean(normalVectors,2))'*normalVectors)); % the tangent plane is constant
    scL = sum(sqrt(sum(diff(ruling3D,1,2).^2))); % the length of the ruling
    sc = (scA+scP)/scL;
    
    % update the first point
    if sc<rulingScore(3,candidateRulings(1,i))
        rulingScore(2,candidateRulings(1,i)) = candidateRulings(2,i);
        rulingScore(3,candidateRulings(1,i)) = sc;
    end
    
    % update the second point -- the score must be minimal
    if sc<rulingScore(3,candidateRulings(2,i))
        rulingScore(2,candidateRulings(2,i)) = candidateRulings(1,i);
        rulingScore(3,candidateRulings(2,i)) = sc;
    end
    
    % update the progress bar
    waitbar(i/n_candidateRulings,hwb);
    set(hTitle,'string',['Ruling Detection Progress: ' num2str(i/n_candidateRulings*100,'%.0f') ' %']);
end
close(hwb)

% filtering the detection : no more than n_maxRulings are allowed for each detection point
indRemove = [];
for i=unique(vect(rulingScore(1:2,:)))'
    [nth,ind] = find(rulingScore(1:2,:)==i);
    if length(ind)>n_maxRulings
        [nth,inds] = sort(rulingScore(3,ind));
        indRemove = [indRemove ; ind(inds(n_maxRulings+1:end))];
    end
end
rulingScore(:,indRemove) = [];

% compute the arc length of the border points
al = [0 cumsum(sqrt(sum(diff(borderPoints,1,2).^2)))];

% convert the ruling indices in arc length
al_rulings = sort(al(rulingScore(1:2,:)));

% uncomment to show the detected rulings
% figure
% plot(corners2D(1,[1:end 1]),corners2D(2,[1:end 1]),'b')
% hold on
% plot(borderPoints(1,rulingScore(1:2,:)),borderPoints(2,rulingScore(1:2,:)),'g')
% title('Detected Rulings')
% end uncomment

%% ruling position curve: exploit the detected rulings to define the score map

n_rulings = size(al_rulings,2);
[mesh_rulings(:,:,1) mesh_rulings(:,:,2)] = meshgrid(linspace(0,1,n_interpolation));
sample_rulings = [vect(mesh_rulings(:,:,1)) vect(mesh_rulings(:,:,2))]';

n_sample_rulings = size(sample_rulings,2);
acc_isolate = zeros(1,n_sample_rulings);
acc_cross = zeros(1,n_sample_rulings);
acc_near = zeros(1,n_sample_rulings);

% to show a progress bar
hwb = waitbar(0,'Cost Function Progress');
hAxes = findobj(hwb,'type','axes');
hTitle = get(hAxes,'title');
cpt=0;
cptmax = n_rulings;

for i=1:n_rulings
    % define the incompatible region for the current ruling
    rightIncompatibleRegion = [al_rulings(1,i) al_rulings(2,[i i]) al_rulings(1,i) ; ...
                               al_rulings(2,[i i]) 1 1];
    leftIncompatibleRegion = [al_rulings(1,i) 0 0 al_rulings(1,i) ; ...
                               al_rulings(2,[i i]) al_rulings(1,[i i])];
    ind_rightIncompatible = (inpolygon(sample_rulings(1,:),sample_rulings(2,:),...
                                 rightIncompatibleRegion(1,:),rightIncompatibleRegion(2,:)));
    ind_leftIncompatible = (inpolygon(sample_rulings(1,:),sample_rulings(2,:),...
                                leftIncompatibleRegion(1,:),leftIncompatibleRegion(2,:)));
    ind_cross = ind_rightIncompatible|ind_leftIncompatible;
    
    % define which rulings are close
	relativeDistance = sqrt(sum((sample_rulings - repmat(al_rulings(:,i),1,n_sample_rulings)).^2));
    ind_near = (relativeDistance<=threshold_closeness);
    ind_leftright = sample_rulings(1,:) > al_rulings(1,i);
    
    % update the counters
    acc_isolate(ind_near) = acc_isolate(ind_near)+1;
    acc_cross(ind_cross&(~ind_near)) = acc_cross(ind_cross&(~ind_near))+1;
    acc_near(ind_near&ind_cross) = acc_near(ind_near&ind_cross) + (2*ind_leftright(ind_near&ind_cross)-1);

    % update the progress bar
     cpt = cpt+1;
     waitbar(cpt/cptmax,hwb);
     set(hTitle,'string',['Cost Function Progress: ' num2str(cpt/cptmax*100,'%.0f') ' %']);
end
close(hwb)


% compute the scores
interpolationRulingScore = reshape((-acc_cross - abs(acc_near))./(acc_isolate),n_interpolation,n_interpolation);
interpolationRulingScoreSign = reshape(sign(acc_near),n_interpolation,n_interpolation);

% uncomment to plot the raw score
% figure
% imagesc(interpolationRulingScore)
% hold on
% grid on
% axis equal
% c = colormap;
% c(1,:) = [1 1 1]; % to avoid the blue background
% colormap(c);
% plot([1 n_interpolation],[1 n_interpolation],'k','linewidth',2)
% nticks = 6;
% set(gca,'Xtick',linspace(1,n_interpolation,nticks),'fontsize',14,...
%         'Xticklabel',num2str(linspace(0,1,nticks)'));
% set(gca,'Ytick',linspace(1,n_interpolation,nticks),...
%         'Yticklabel',num2str(linspace(0,1,nticks)'));
% set(gca,'Xlim',[1 n_interpolation]);
% set(gca,'Ylim',[1 n_interpolation]);
% axis xy
% xlabel('s1')
% ylabel('s2')
% end uncomment


% extract the ruling position curve from the score map
% preprocessing the score map
interpolationRulingScore(triu(true(n_interpolation))) = nan;
interpolationRulingScoreSign(triu(true(n_interpolation))) = nan;

% extract the curves
ind = ~(isnan(interpolationRulingScore)|isinf(interpolationRulingScore));
interpolationRulingScore = interpolationRulingScore - min(interpolationRulingScore(ind));
interpolationRulingScore = interpolationRulingScore/max(interpolationRulingScore(ind));
interpolationRulingScore(~ind) = 0;

% apply the threshold to define the regions
th_region = graythresh(interpolationRulingScore(ind)); 
mask_region = imfill(double(interpolationRulingScore>th_region));
interpolationRulingScoreSign = interpolationRulingScoreSign.*mask_region;

% thin the mask to avoid contour detection
mask_region = imerode(mask_region,strel('disk',round(0.025*n_interpolation)));
ind = (interpolationRulingScoreSign==0)&(mask_region==1);
opt = 'nearest';
[iref,jref] = find((interpolationRulingScoreSign~=0)&(mask_region==1));
[i,j] = find(ind);
interpolationRulingScoreSign(ind) = griddata(iref,jref,interpolationRulingScoreSign((interpolationRulingScoreSign~=0)&(mask_region==1)),i,j,opt);
rulingPositionCurves = edge(uint8(interpolationRulingScoreSign),'sobel');

% keep only the significant part of the mask
mask_region_bw = bwlabel(mask_region==1);
rulingsLabel = mask_region_bw(sub2ind(size(mask_region_bw),round((n_interpolation-1)*al_rulings(2,:)+1),round((n_interpolation-1)*al_rulings(1,:)+1)));
newmask = zeros(size(mask_region));
for i=1:max(mask_region_bw(:))
    ind = rulingsLabel==i;
    newmask(round(min((n_interpolation-1)*al_rulings(2,ind)+1)):round(max((n_interpolation-1)*al_rulings(2,ind)+1)),...
            round(min((n_interpolation-1)*al_rulings(1,ind)+1)):round(max((n_interpolation-1)*al_rulings(1,ind)+1))) = 1;
end
mask_region = newmask.*mask_region_bw~=0;

% uncomment to show the detected interesting regions
% figure
% imagesc(mask_region)
% hold on
% grid on
% axis equal
% colormap([1 1 1 ; 0 0 0]);
% plot([1 n_interpolation],[1 n_interpolation],'k','linewidth',2)
% set(gca,'Xtick',linspace(1,n_interpolation,nticks),'fontsize',14,...
%         'Xticklabel',num2str(linspace(0,1,nticks)'));
% set(gca,'Ytick',linspace(1,n_interpolation,nticks),...
%         'Yticklabel',num2str(linspace(0,1,nticks)'));
% set(gca,'Xlim',[1 n_interpolation]);
% set(gca,'Ylim',[1 n_interpolation]);
% axis xy
% xlabel('s1')
% ylabel('s2')
% title('detected regions')
% end uncomment


% extract the ruling position curve
interpolationRulingScore(isnan(interpolationRulingScore)) = 0;
interpolationRulingScoreSign(isnan(interpolationRulingScoreSign)) = 0;
rulingPositionCurves = rulingPositionCurves.*mask_region; % to eliminate the region borders (keep only 


% plot the detected curve 1/2
% figure
% hold on
% grid on
% axis equal
% plot([1 n_interpolation],[1 n_interpolation],'k','linewidth',2)
% set(gca,'Xtick',linspace(1,n_interpolation,nticks),'fontsize',14,...
%         'Xticklabel',num2str(linspace(0,1,nticks)'));
% set(gca,'Ytick',linspace(1,n_interpolation,nticks),...
%         'Yticklabel',num2str(linspace(0,1,nticks)'));
% set(gca,'Xlim',[1 n_interpolation]);
% set(gca,'Ylim',[1 n_interpolation]);
% axis xy
% xlabel('s1')
% ylabel('s2')
% uncomment the next one as well

% enforce monotonicity of each curve
labeledRulingPositionCurves = bwlabel(rulingPositionCurves);
n_curves = max(labeledRulingPositionCurves(:));
labeledRulingPositionCurves_old = labeledRulingPositionCurves;
labeledRulingPositionCurves = zeros(n_interpolation);
for i=1:n_curves
    [x,y] = find(labeledRulingPositionCurves_old == i);
    xyc = unique([-y x],'rows');
    x = xyc(:,2);
    y = -xyc(:,1);
    
    % plot the detected curve 2/2
%     plot(y,x,'k','linewidth',3)
    % end uncomment
    
    [x,y] = find(labeledRulingPositionCurves_old == i);
    [x_new,y_new] = enforceDecreasing(x,y);
    if length(x_new>=n_rulingCurveMin)
        indRMXY = round(x_new)<1 | round(x_new)>n_interpolation | round(y_new)<1 | round(y_new)>n_interpolation;
        x_new(indRMXY) = [];
        y_new(indRMXY) = [];
        labeledRulingPositionCurves(sub2ind([n_interpolation n_interpolation],round(x_new),round(y_new)))=1;
    end
end
labeledRulingPositionCurves = bwlabel(labeledRulingPositionCurves);

% enforce global compatibility and guess the ruling positions where the detection failed
labeledRulingPositionCurves = cutAndGuess(labeledRulingPositionCurves); % to be improved

% filter the ruling position curve
n_f = 5;
f = ABFilters1D('gauss_x',3*n_interpolation/512,n_f);
labeledRulingPositionCurves_new = zeros(size(labeledRulingPositionCurves));
for i=1:max(labeledRulingPositionCurves(:))
    [ic,jc] = find(labeledRulingPositionCurves==i);
    ijc = unique([ic -jc],'rows');
    ic = ijc(:,1);
    jc = -ijc(:,2);

    % apply the filter
    ic = conv([repmat(ic(1),n_f,1) ; ic ; repmat(ic(end),n_f,1)],f);
    ic = ic(2*n_f+1:end-2*n_f);
    jc = conv([repmat(jc(1),n_f,1) ; jc ; repmat(jc(end),n_f,1)],f);
    jc = jc(2*n_f+1:end-2*n_f);
    
    labeledRulingPositionCurves_new(sub2ind(size(labeledRulingPositionCurves_new),round(ic),round(jc))) = i;
    
end
labeledRulingPositionCurves = labeledRulingPositionCurves_new;
% enforce global compatibility and guess the ruling positions where the detection failed
labeledRulingPositionCurves = cutAndGuess(labeledRulingPositionCurves); % to be improved


% uncomment to show the ruling position curve
% figure
% axis equal
% hold on
% plot([1 n_interpolation],[1 n_interpolation],'k','linewidth',2)
r = [];
c = [];
for i=1:max(labeledRulingPositionCurves(:))
    [ic,jc] = find(labeledRulingPositionCurves==i);
    ijc = unique([ic -jc],'rows');
    ic = ijc(:,1);
    jc = -ijc(:,2);
    
    ic = conv([repmat(ic(1),n_f,1) ; ic ; repmat(ic(end),n_f,1)],f);
    ic = ic(2*n_f+1:end-2*n_f);
    jc = conv([repmat(jc(1),n_f,1) ; jc ; repmat(jc(end),n_f,1)],f);
    jc = jc(2*n_f+1:end-2*n_f);
    r = [r ; ic];
    c = [c ; jc];
    
    % uncomment to show
%     plot(jc,ic,'k','linewidth',3)
end
% uncomment to show
% axis equal
% grid on
% plot([1 n_interpolation],[1 n_interpolation],'k','linewidth',2)
% set(gca,'Xtick',linspace(1,n_interpolation,nticks),'fontsize',14,...
%         'Xticklabel',num2str(linspace(0,1,nticks)'));
% set(gca,'Ytick',linspace(1,n_interpolation,nticks),...
%         'Yticklabel',num2str(linspace(0,1,nticks)'));
% set(gca,'Xlim',[1 n_interpolation]);
% set(gca,'Ylim',[1 n_interpolation]);
% axis xy
% xlabel('s1')
% ylabel('s2')
% grid on
% axis equal
% cm = colormap;
% cm(1,:) = [1 1 1]; % to avoid the blue background
% colormap(cm);
% plot([1 n_interpolation],[1 n_interpolation],'k','linewidth',2)
% set(gca,'Xtick',linspace(1,n_interpolation,nticks),'fontsize',14,...
%         'Xticklabel',num2str(linspace(0,1,nticks)'));
% set(gca,'Ytick',linspace(1,n_interpolation,nticks),...
%         'Yticklabel',num2str(linspace(0,1,nticks)'));
% set(gca,'Xlim',[1 n_interpolation]);
% set(gca,'Ylim',[1 n_interpolation]);
% axis xy
% xlabel('s1')
% ylabel('s2')
% end uncomment

%%
%% Paper initialization

% flat paper
pp.r = dx/dy;
alpp = unique([r c],'rows')';

% pre-processing to avoid crossing
for k=alpp(:)'
    [i,j] = find(alpp==k);
    if length(i)>1 % if there is more than one point
        % a different traitment is recquired according to the relative position of the points (first or
        % second point of the ruling
        % first points
        i1 = i(i==1);
        j1 = j(i==1);
        lp = alpp(sub2ind(size(alpp),3-i1,j1)); % find the points
        [nth,ind] = sort(lp,'descend'); 
        j1 = j1(ind); % sort the points
        i1 = i1(ind);
        stp = linspace(-.5,0,length(j1)+2);
        alpp(sub2ind(size(alpp),i1,j1)) = k + stp(2:end-1); % spread the points
        
        % second points
        i2 = i(i==2);
        j2 = j(i==2);
        lp = alpp(sub2ind(size(alpp),3-i2,j2)); % find the points
        [nth,ind] = sort(lp,'descend'); 
        j2 = j2(ind); % sort the points
        i2 = i2(ind);
        stp = linspace(0,.5,length(j2)+2);
        alpp(sub2ind(size(alpp),i2,j2)) = k + stp(2:end-1); % spread the points
    end
end

pp.al = alpp/(n_interpolation+1);
pp.be = zeros(1,size(pp.al,2));
% remove ruling intersections
[pp.al , pp.be] = uncross(pp.al , pp.be);

paper = newPaperFast(pp);

% angle computation
pp = paper.ppnew;
pp = paperComputeBe(pp,paper.sh,points2D_all,points3D_all);
pp.al(:,pp.be==0) =[];
pp.be(:,pp.be==0) =[]; % something wrong probably happened

%%
%% guiding ruling selection by curve fitting

paperSubOpt = newPaperFast(pp);
ppNew = paperOptimalParameterisationWB(paperSubOpt,th_curveFitting);
ppNew = paperComputeBe(ppNew,paper.sh,points2D_all,points3D_all);

% uncomment to show the flat paper
figure(33);
plotFlatPaper(paperSubOpt,3);
axis equal
axis off
% end uncomment

paper = newPaperFast(ppNew);

% uncomment to show the paper with the selected guiding rulings
plotFlatPaper(paper);
axis equal
% end uncomment

% rigid transformation
[m(:,:,1),m(:,:,2)] = meshgrid(px,py);
Pm = reshape(m,nx*ny,2);
W = GTPSW_TPS_kernelize(points2D_all',Pm);
PM = TPS*W';
M = reshape(PM',[ny,nx,3]);
[Mpaper,rigidTransformation] = meshRigidTransformation(paperMesh(paper,nx,ny),M);
MpaperSubOpt = meshRigidTransformation(paperMesh(paperSubOpt,nx,ny),M);

interp.origin = MpaperSubOpt;
interp.destination = Mpaper;

% uncomment to show the initial surface, the paper before the guiding ruling selection and the paper
% after the guiding ruling selection.
figure(20);
plotMesh(M,'k')
hold on
plotMesh(MpaperSubOpt,'b')
plotMesh(Mpaper,'g')
axis equal
plot3(points3D_all(1,:),points3D_all(2,:),points3D_all(3,:),'b.','markersize',10)
title('initialisation of the paper');
legend({'Initial surface', 'Paper:pre-guiding-ruling-select', 'Paper:post-guiding-ruling-select'}, 'Interpreter', 'Latex', 'Location','southoutside');
pause(0.1);
% end uncomment

