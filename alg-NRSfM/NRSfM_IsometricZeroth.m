function [resultNRSfM, numEdges] = NRSfM_IsometricZeroth(Data, nng, K, varargin)

    p = inputParser;
    expectedDeformationModel = {'mdh-exactly-isometric', 'mdh-extensible-isometric', 'mdh-extensible-upto-equiareality', 'mdh-inextensible-isometric'};
    expectedSolvers = {'mosek', 'sdpt3', 'yalmip', 'sedumi'};
    defaultDeformationModel = expectedDeformationModel{2};
    defaultLambda = 1;
    defaultEpsilon = 0.00;
    defaultSolver = expectedSolvers{1};
    defaultScaleCompensation = true;

    addRequired(p,'Data');
    addRequired(p,'nng');
    addRequired(p,'K');
    addParameter(p,'DeformationModel',defaultDeformationModel,@(x) any(validatestring(x,expectedDeformationModel)));
    addParameter(p,'AreaTensorTraceLambda',defaultLambda);
    addParameter(p,'AreaMinimisationLambda',defaultLambda);
    addParameter(p,'Solver',defaultSolver,@(x) any(validatestring(x,expectedSolvers)));
    addParameter(p,'EPSILON',defaultEpsilon);
    addParameter(p,'RescaleData',defaultScaleCompensation);

    parse(p,Data, nng, K, varargin{:});
    
    trI = TriangleIndexer(nng);
    edI = EdgeIndexer(nng);
    dtI = DeltaTensorIndexer(nng);

    params = p.Results;

    % for isometric:
    % params.IsometryCostLambda = 40000000;
    % params.DeltaTraceLambda = 0.001;
    % params.MdhLambda = 90;

    params.IsometryCostLambda = 40000000;
    params.DeltaTraceLambda = 0.001;
    params.MdhLambda = 1000;


    [nFiles, nPts, kN] = checkDataSanity(params);

    N = size(Data.p,2);
    M = size(Data.p(1).p,2);
    Kinv = inv(K);  
    for iImg = 1:N
        directionVectors = [];
        for iPts = 1:M
            dirVec_j = Kinv*[Data.p(iImg).p(1:2, iPts); 1];
            dirVec_j = dirVec_j./norm(dirVec_j);
            directionVectors = [directionVectors dirVec_j];
        end
        directionVector{iImg} = directionVectors;
    end
    
    
    cvx_begin
    cvx_precision low
    if(strcmp(params.Solver,expectedSolvers{1}))
        cvx_solver mosek
    elseif(strcmp(params.Solver,expectedSolvers{2}))
        cvx_solver sdpt3
    elseif(strcmp(params.Solver,expectedSolvers{3}))
        cvx_solver yalmip
    elseif(strcmp(params.Solver,expectedSolvers{4}))
        cvx_solver sedumi
    end

   
    variable geod(edI.numberOfEdges) nonnegative;

    if(strcmp(params.DeformationModel,expectedDeformationModel{3}))
         variable deltaD(nFiles,(nPts+dtI.numberOfEdges), (nPts+dtI.numberOfEdges)) nonnegative;
        variable A(trI.numberOfTriangles) nonnegative;
        variable KK;
        varSize = numel(deltaD) + numel(geod) + numel(A);
        printing_frequency = 0.6;
        areaObj = 0;
    else
         variable deltaD(nFiles,nPts, nPts) nonnegative;
        varSize = numel(deltaD) + numel(geod);
        A = zeros(trI.numberOfTriangles);
        X = zeros(nFiles, dtI.numberOfEdges, dtI.numberOfEdges);
        printing_frequency = 0.75;
    end
    
    numEdges = dtI.numberOfEdges;


    mdhObjective = 0;
    isometryObjective = 0;
    deltaTraceObjective = 0;
    areaTensorTraceObjective = 0;

    numTotTriangles = 0;
    numTotEdges = 0;

    subject to:
    for iF = 1:nFiles
        for iP = 1:nPts
            vecJ = directionVector{iF}(:,iP);
            mdhObjective = mdhObjective + pow_p(deltaD(iF, iP, iP), -1);
            for iN = 1:(kN)

                nN1 = nng(iP, iN);
                edgeIndex = edI.getIndex(iP, nN1);
                
                vecQ = directionVector{iF}(:,nN1);

                dVec = deltaD(iF, iP, iP) + deltaD(iF, nN1, nN1) - (2*deltaD(iF, iP, nN1)*dot(vecJ, vecQ));

                if ((iN < kN) && (strcmp(params.DeformationModel,expectedDeformationModel{3})))
                    nN2 = nng(iP, iN+1);
                    trIndex = trI.getIndex(iP, nN1, nN2);
                    
                    vecR = directionVector{iF}(:,nN2);
                    
                    e_jq = nPts + dtI.getIndex(iP, nN1);
                    e_qr = nPts + dtI.getIndex(nN1, nN2);
                    e_jr = nPts + dtI.getIndex(iP, nN2);

                    x = [deltaD(iF, e_jq, e_jq) deltaD(iF,e_jq, e_jr) deltaD(iF, e_jr, e_jr) deltaD(iF,e_jq,e_qr) deltaD(iF, e_qr, e_jr) deltaD(iF, e_qr, e_qr)];
                    [G] = getQuadraticCoefficientsTriangleArea(vecJ, vecQ, vecR);
                    areaObj = areaObj + abs((1/4)*G*x.' - A(trIndex));
                    

                    deltaD(iF, iP, nN1)^2 <= x(1); 
                    deltaD(iF, iP, nN2)^2 <= x(3); 
                    deltaD(iF, nN1, nN2)^2 <= x(6); 

                    numTotTriangles = numTotTriangles + 1;
                end

                if ( (strcmp(params.DeformationModel,expectedDeformationModel{3})) || (strcmp(params.DeformationModel,expectedDeformationModel{2})) )
                    isometryObjective = isometryObjective + abs(dVec - geod(edgeIndex) );
                elseif(strcmp(params.DeformationModel,expectedDeformationModel{1}))
                    dVec == geod(edgeIndex); 
                else
                    dVec <= geod(edgeIndex); 
                end

                numTotEdges = numTotEdges + 1;

            end
        end

        B = squeeze(deltaD(iF,:,:));
        n = size(B,1);
        B == semidefinite(n); 
        deltaTraceObjective = deltaTraceObjective + trace(B);
    end

    
    if(strcmp(params.DeformationModel,expectedDeformationModel{3}))
       
        sum(sum(geod)) == 1;
        A >= params.EPSILON;
        geod >= params.EPSILON;
        
        sum(sum(A)) <= sum(sum(geod));
    else
        geod >= params.EPSILON; 
        sum(sum(geod)) == 1; 
    end

    mdhObjective = mdhObjective/(nFiles*nPts);
    areaTensorTraceObjective = areaTensorTraceObjective/numTotTriangles;
    isometryObjective = isometryObjective/numTotEdges;
    deltaTraceObjective = deltaTraceObjective/(nFiles*nPts);
    

    if(strcmp(params.DeformationModel,expectedDeformationModel{1}) || strcmp(params.DeformationModel,expectedDeformationModel{4}))
        minimize( (params.MdhLambda*mdhObjective) + (params.DeltaTraceLambda*deltaTraceObjective) );
    elseif(strcmp(params.DeformationModel,expectedDeformationModel{2}))
        minimize( (params.MdhLambda*mdhObjective) + (params.IsometryCostLambda*isometryObjective)  + (params.DeltaTraceLambda*deltaTraceObjective));
    else
        minimize( (params.MdhLambda*mdhObjective) + (params.IsometryCostLambda*isometryObjective)  + (params.DeltaTraceLambda*deltaTraceObjective) ...
           +  (params.AreaMinimisationLambda*areaObj) ); 
    end
    cvx_end


    depthOfPoints = zeros(nFiles, nPts);

    for iF = 1:nFiles
        gramMatrixList(iF).G = full(squeeze(deltaD(iF,:,:)));

        [U, S, V] = svd(squeeze(deltaD(iF,:,:)));
        D = U(:,1)*S(1,1)*V(:,1).';
        depthOfPoints(iF,:) = sqrt(diag(D(1:nPts,1:nPts)));
    end


    if(strcmp(params.DeformationModel,expectedDeformationModel{3}))
        areaMatrix = A;
    end
    geodesicMatrix = geod;

    
    if(params.RescaleData == true)
        [reconstruction, scaleList] = getReconstruction(Data, depthOfPoints, directionVector, nFiles);
    else
        [reconstruction, scaleList] = getReconstructionUnscaled(Data, depthOfPoints, directionVector, nFiles);
    end
     

    resultNRSfM.status = cvx_status;
    resultNRSfM.reconstruction = reconstruction;
    resultNRSfM.scaleList = scaleList;
    resultNRSfM.geodesicMatrix = geodesicMatrix;
    if(strcmp(params.DeformationModel,expectedDeformationModel{3}))
        resultNRSfM.areaMatrix = areaMatrix;
    end
    resultNRSfM.depthOfPoints = depthOfPoints;
    resultNRSfM.summedGeodSquared = sum(geod.^2);
    resultNRSfM.gramMatrices = gramMatrixList;
    

end

function [nFiles, nPts, kN] = checkDataSanity(params)
    nFiles = size(params.Data.Pgth,2);
    nPts = size(params.Data.Pgth(1).P,2);
    kN = size(params.nng,2);    
end
