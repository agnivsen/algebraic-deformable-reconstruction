function [Data, nng, K, areaTemplate, flatTemplate] = syntheticDataGenerator(nFiles, gridSize, nNeighbors, noise, showDebugPlots, equiarealityStrength)
   %% DataEquiarealNRSfMFormatter creates randomly deforming equiareal surfaces based on the code from [Perriollat et. al., 2013]
    % Following the paper model of [Perriollat et. al., 2013], the generated _isometric_ surface is randomly perturbed and then brought
    % back to an equiareal configuration using a non-linear least squares (_Levenberg-Marquardt_)
    %
    % *PARAMETER*:
    %
    % *  _nFiles_: number of files/images required (n) [minimum = 2]
    % * _gridSize_: [m = m1 x m2] grid of correspondences required, where (m1, m2 >=1)
    % * _nNeighbors_: number of neighbors (kN) [minimum = 2, maximum = (m-1)]
    % * _noise_: scaling factor for additive noise from uniform random distribution, to be added on 2D point correspondences on the image plane,
    %                   in pixel units
    % * _showDebugPlots_: boolean indicating if displaying debug plots are required
    % * _equiarealityStrength_: a variable in [0,1], 0 indicating completely _isometric_ while 1 indicates compteley _equiareal_
    %
    % *RETURNS*
    %
    % * _Data_: A *struct* containing two subfields Data.p and Data.Pgth,
    %                   where Data.p(i).p in R^[3 x m] = homogeneous pixel coordinates
    %                   and    Data.Pgth(i).P in R^[3 x m] = 3D GT coordinates
    % * _nng_: [m x kN] array of neighbors
    % * _K_: [3 x 3] camera intrinsics, randomly generated
    % * _areaTemplate_: [m x (kN-1)] area of triangles on the flattened template, fromed by vertices [P_j, P_q, P_(q+1)] for all j in [1,m] and q in [1,(NNG(j)-1)]
    %
    %  [Perriollat et. al., 2013]: Perriollat, Mathieu, and Adrien Bartoli. "A computational model of bounded developable surfaces with application to image based three dimensional reconstruction." Computer Animation and Virtual Worlds 24.5 (2013): 459-476.
    %                                             Code at: http://igt.ip.uca.fr/~ab/code_and_datasets/index.php

    close all;
    
    pathToStDRv1p0c = './Utils/St-DR_v1p0c/'; 
    
    modelPath = strcat(pathToStDRv1p0c, 'model/');
    miscPath =  strcat(pathToStDRv1p0c, 'misc/');
    otherPath =  strcat(pathToStDRv1p0c, 'other/');
    
    addpath(pathToStDRv1p0c);
    addpath(modelPath);
    addpath(miscPath);
    addpath(otherPath);

    if~exist('equiarealityStrength','var')
        equiarealityStrength = 0.1;
    end
    
    
    assert(numel(gridSize) == 2,'DataEquiarealNRSfMFormatter: for <gridSize>, expecting a [1 x 2] array of dimensions');
    assert(isequal(size(gridSize,1),1),'DataEquiarealNRSfMFormatter: for <gridSize>, expecting a [1 x 2] array of dimensions');
    assert(gridSize(1) >= 1 && isinteger(int8(gridSize(1))),'DataEquiarealNRSfMFormatter: <gridSize> dimension should be integer and >= 1');
    assert(gridSize(2) >= 1 && isinteger(int8(gridSize(2))),'DataEquiarealNRSfMFormatter: <gridSize> dimension should be integer and >= 1');
    
    assert(nFiles > 1 && isinteger(int8(nFiles)), 'DataEquiarealNRSfMFormatter: <nFiles> should be an integer and greater than 1');

%     assert((equiarealityStrength >=0 ) && (equiarealityStrength <= 1),... 
%                                         'DataEquiarealNRSfMFormatter: controls how much the data differs from isometric model, 0 is actually isometric. Needs to be in range [0,1]')
    
    
    if exist('nNeighbors','var')
        assert(nNeighbors>1,'DataEquiarealNRSfMFormatter: at least 2 neighbors needed for equiareal data generation using this method');
        assert(nNeighbors<(gridSize(1)*gridSize(2)),'DataEquiarealNRSfMFormatter: no. of neighbors must be less than total no. of pt.s in every image');
        nK = nNeighbors;
    else
        nK = 3;
    end

    showDbg = false;
    if exist('showDebugPlots','var')
        if(showDebugPlots == true)
            showDbg = true;
        end
    end
    
    if( ~(exist('noise','var')) || (numel(noise) == 0))
        noise = 0;
    else
        fprintf('NOISE added to <strong>correspondences</strong> = %d pixels [multiplier to uniform distribution]\n\n', noise);
    end
    
%     fx = round(random_number_within_range(700,900,1));
%     fy = round(random_number_within_range(700,900,1));
    
    fx = 800;
    fy = 800;
    
    IMG_HEIGHT = 1080*2;
    IMG_WIDTH = 1921*2;
    
    cx = IMG_WIDTH/2;
    cy = IMG_HEIGHT/2;
    
    K = [fx 0 cx; 0 fy cy; 0 0 1];
    
    
    % generate a paper parameter structure
    pp1.type = 'al_be'; % the id of the parameterization
    pp1.r = 1.4;%random_number_within_range(1,2,1); % width length ratio
%     arc1 = random_number_within_range(0, 1, 3);
%     arc2 = random_number_within_range(1, 2, 3);
%     pp.al = [arc1; arc2];
    pp1.al = [0.57 0.6374 0.9118 ; 1.1469 0.8737 1.111]; % arc lengths of the rulings
    
    
    pp1.be = [0 0 0]; % normalized bending angles

    % interpolate the rulings
    pp1 = paperParameterisationConversion(pp1,'al_be',5);

    % compute the paper
    paper = newPaperFast(pp1);

%     % display the flat paper
%     plotFlatPaper(paper);

    % evaluate the 3D mesh
    templatePaper = paperMesh(paper,gridSize(1),gridSize(2));
    
    [template3D] = meshToPointcloud(templatePaper);
    
    flatTemplate = template3D;
    
    tempdata.p(1).p = template3D.';
        
    [nng] = getNeighborhoodDuplicated(tempdata,nK);

    [areaTemplate] = computePiecewiseAreaOfMesh(template3D, nng);
    
    if(showDbg)
        figure;
    end
    
    for ii = 1:size(nng,1)
        for jj = 1:nK
            nn = nng(ii,jj);
            if(showDbg)
                if(jj==1)
                    plot([template3D(ii,1) template3D(nn,1)].',[template3D(ii,2) template3D(nn,2)].','r--','LineWidth',2); hold on;
                else
                    plot([template3D(ii,1) template3D(nn,1)].',[template3D(ii,2) template3D(nn,2)].','LineWidth',2); hold on;
                end
            end
            dist = sqrt( (template3D(ii,1) - template3D(nn,1)).^2 + (template3D(ii,2) - template3D(nn,2)).^2 );
            geodesic{ii,nn} = dist;
        end
    end
    
    
        % generate a paper parameter structure
    pp.type = 'al_be'; % the id of the parameterization
    pp.r = 1.4;%random_number_within_range(1,2,1); % width length ratio
%     arc1 = random_number_within_range(0, 1, 3);
%     arc2 = random_number_within_range(1, 2, 3);
%     pp.al = [arc1; arc2];
    pp.al = [0.57 0.6374 0.9118 ; 1.1469 0.8737 1.111]; % arc lengths of the rulings
   
    
    for iFile = 1:nFiles
    
        bend = random_number_within_range(-15, 15, 3);
        pp.be = bend; % normalized bending angles

        % interpolate the rulings
        ppn = paperParameterisationConversion(pp,'al_be',5);

        % compute the paper
        paper = newPaperFast(ppn);

%         % display the flat paper
%         plotFlatPaper(paper);

        % evaluate the 3D mesh
        mesh3D = paperMesh(paper,gridSize(1),gridSize(2));
        
        [pointcloud] = meshToPointcloud(mesh3D);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Uncomment the next line of code (and comment the previous one)
        %%% to generate 1D data that lies on a straightline
        
%         [pointcloud] = getStraightLine(gridSize(1)*gridSize(2));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        C = mean(pointcloud);
        
        pointcloud = pointcloud.' - C.';
        pointcloud = pointcloud.';
%         pointcloud = pointcloud + 0.01.*rand(size(pointcloud));
        
        rot1 = random_number_within_range(-10,10,1);
        rot2 = random_number_within_range(-10,10,1);
        rot3 = random_number_within_range(-45,45,1);
        
        rX = rotx(deg2rad(rot1)); rY = roty(deg2rad(rot2)); rZ = rotz(deg2rad(rot3));
        
        Rs = rX*rY*rZ;
        
        tx = random_number_within_range(-0.2,0.2,1);
        ty = random_number_within_range(-0.2,0.2,1);
        tz = random_number_within_range(0.25,0.35,1);
        
        t = [tx; ty; tz];
        
        pointcloud = Rs*pointcloud.' + t;
        
        pointcloud = pointcloud.';    

        isometricCloud = pointcloud;

        if(showDbg)
            fPC = figure;
            subplot(2,1,1);
            set(gcf,'Color','w'); hold on;
            set(gca,'Color','w'); hold on;
            pcshow(pointcloud,'r','MarkerSize', 40); hold on;
        end

        low = min(pointcloud);
        high = max(pointcloud);
        lengthHighest = EuclideanDistance(low, high);
        multiplier = equiarealityStrength * lengthHighest / 0.5;

        pointcloud = pointcloud + multiplier.*(rand(size(pointcloud)) - 0.5);

        randomCloud = pointcloud;

%         pcshow(pointcloud, 'g', 'MarkerSize', 30); hold on;

        if(equiarealityStrength > 0.01)

            [X] = serializePoints(pointcloud);

            [E, areaErrorInitial] = areaError(X, areaTemplate, nng, randomCloud);

            tolerance = 10.^(-17);
            options = optimoptions('lsqnonlin','Display','off', 'StepTolerance', tolerance,'Algorithm','levenberg-marquardt', ...
                                                            'FunctionTolerance', tolerance, 'OptimalityTolerance', tolerance, 'MaxIterations', 150, ...
                                                            'MaxFunctionEvaluations', 10000);   

            [X] = lsqnonlin(@areaError,X,[], [], options, areaTemplate, nng, randomCloud);

            [pointcloud] = deserializePoints(X);
        

            if(showDbg)
                pcshow(pointcloud, 'b', 'MarkerSize', 40); hold on;

                legend('Original isometric',  'Equiareal - refined');

                [E, areaErrorFinal] = areaError(X, areaTemplate, nng, randomCloud);
                subplot(2,1,2);
                set(gcf,'Color','w');
                plot(areaErrorFinal,'k-','LineWidth', 2); hold on;

                pause(1);
                close(fPC);
            end
        end
        
        imgPoints = [];
        
        noiseProfile = [];
        for ii=1:size(pointcloud,1)
            P_ = pointcloud(ii,:);
            n1 = noise*rand;
            px = P_(1)/P_(3)*fx + cx + n1;
            n2 = noise*rand;
            py = P_(2)/P_(3)*fy + cy + noise*rand;
            noiseProfile = [noiseProfile n1 n2];
            p = [px py 1];
            imgPoints = [imgPoints;p];
        end
        
        if(abs(noise) > 10^(-3))
            fprintf('====> Noise profile for img. %d: [<strong>mean</strong> = %d], [<strong>median</strong> = %d], [<strong>SD</strong> = %d]\n',iFile, mean(noiseProfile), median(noiseProfile),std(noiseProfile));
        end
        
        if(showDbg)
            subplot(1,2,1);
            scatter(imgPoints(:,1), imgPoints(:,2), 'r*');
            xlim([0 IMG_WIDTH]);
            ylim([0 IMG_HEIGHT]);
            titStr = strcat('Frame no.', num2str(iFile));
            title(titStr);
        end

        
        Data.p(iFile).p = imgPoints.';
        Data.Pgth(iFile).P = pointcloud.';
        Data.v = ones(nFiles, (gridSize(1) * gridSize(2)));
        
        if(showDbg)
            subplot(1,2,2);
            pcshow(pointcloud,'MarkerSize', 70); hold on;
            set(gcf,'color','w');
            set(gca,'color','w');
            pause(1);
        end
    
    end

% % % Following code block makes points invisible    
% % %     toHideCount = 10;
% % %     start = 1;
% % %     
% % %     for iF = 1:nFiles
% % %         fprintf('<strong>Hidden:</strong> for img. no. %d\n', iF)
% % % %         toHide = randperm((gridSize(1) * gridSize(2)));
% % %         toHide = linspace(1,(gridSize(1) * gridSize(2)), (gridSize(1) * gridSize(2)));
% % %         toHide = toHide(start:(start+toHideCount-1));
% % %         if((start+toHideCount) < (gridSize(1) * gridSize(2)))
% % %             start = start + 1;
% % %         end
% % %         for ii = 1: toHideCount
% % %             fprintf('%d  ', toHide(ii));
% % %         end
% % %         Data.v(iF,toHide) = 0;
% % %         fprintf('\n');    
% % %     end
    
    
    if(showDbg)
        close all;
    end

end

function [E, areaError] = areaError(X, areaTemplate, nng, randomCloud)
    [pointcloud] = deserializePoints(X);
    area = computePiecewiseAreaOfMesh(pointcloud, nng);
    E = [];

    W1 = 0.7;
    W2 = 0.1;

    pushBackward = 0;

    for ii = 1:size(area,1)
        for  jj = 1:size(area,2)
            error = abs(area{ii,jj}^2 - areaTemplate{ii,jj}^2);
            E = [E W1*error];
            pushBackward = pushBackward + (1/area{ii,jj});
        end
    end

    areaError = E;

%     [smoothness] = computeLaplacianSmoothnessMesh(pointcloud, nng);
%     E = [E W2.*smoothness];

    rms = RMSE3D(pointcloud, randomCloud);
    E = [E (1 - W1)*rms];
end

function [X] = serializePoints(pointcloud)
    nP = size(pointcloud,1);
    X = reshape(pointcloud,[3*nP,1]);
end

function [pointcloud] = deserializePoints(X)
    nP = numel(X)/3;
    pointcloud = reshape(X,[nP,3]);
end


function [pointcloud] = getStraightLine(numPts)
    x = linspace(-1, 1, numPts);
    pointcloud = zeros(numPts, 3);
    pointcloud(:,1) = x.';
end

function [pointcloud] = meshToPointcloud(mesh3D)

    numPts = size(mesh3D,1)*size(mesh3D,2);
    
    pointcloud = zeros(numPts, 3);
    
    index = 1;
    
    for ii = 1: size(mesh3D,1)
        for jj = 1:size(mesh3D,2)
            pointcloud(index,:) = [mesh3D(ii,jj, 1) mesh3D(ii,jj, 2) mesh3D(ii,jj, 3)];
            index = index + 1;
        end
    end
    
    
end
