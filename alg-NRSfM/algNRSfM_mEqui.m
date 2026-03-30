function [reconstruction, solutionStatus] = algNRSfM_mEqui(Data, K, nng, odReconstruction, sumSqGeod, eps, solver, lambda)

    if ~exist('eps', 'var')
        EPSILON = 0.1;
    else
        EPSILON = eps;
    end
    
    if ~exist('solver', 'var')
        solver = 'sdpt3';
    end

    if ~exist('lambda', 'var')
        lambda = 10;
    end

    N = size(Data.p,2);
    M = size(Data.p(1).p,2);
    Kinv = inv(K);  
    edI = EdgeIndexer(nng);
    depthOfPointsOD = zeros(N,M);


    for iImg = 1:N
        directionVectors = []; depths = [];
        for iPts = 1:M
            dirVec_j = Kinv*[Data.p(iImg).p(1:2, iPts); 1];
            dirVec_j = dirVec_j./norm(dirVec_j);
            directionVectors = [directionVectors dirVec_j];
            depths = [depths norm(odReconstruction{iImg}(iPts,:))];
        end
        depthOfPointsOD(iImg,:) = depths + EPSILON;
        directionVectorsList{iImg} = directionVectors;
    end

    % assert(max(max(depthOfPointsOD)) < 1, 'Max. depth from OD prior needs to be less than 1');
    % If OD prior does have depths > 1, we need to rescale the OD prior appropriately

    trI = TriangleIndexer(nng);

    assert(trI.numberOfTriangles > 1, 'There has to be at least 1 valid triangle in the NNG to impose equiareal constraint!');

    fprintf(['\n(Equiareal NRSfM has <strong>%d triangles</strong> with <strong>%d total</strong> edges, <strong>%d utilised</strong> edges and <strong>moment matrix</strong> of size <strong>[%d x %d]</strong>' ...
        ' in this data).\n\n'], trI.numberOfTriangles, trI.numberOfEdges, trI.numUsedEdges, M+1+trI.numUsedEdges, M+1+trI.numUsedEdges);
    
    lambda_equi = lambda/10000;
    lambda = lambda/10;
    edgeVisits = zeros(N, trI.numUsedEdges);

    for iImg = 1:N
        for ii = 1:M
            pInfty_p = norm(odReconstruction{iImg}(ii,:)) + EPSILON;
    
            for iNeighbor = 1:size(nng, 2)
                iN = nng(ii, iNeighbor);
                pInfty_q = norm(odReconstruction{iImg}(iN,:)) + EPSILON;
                scalarProd = dot(directionVectorsList{iImg}(:,ii), directionVectorsList{iImg}(:,iN));
    
                constTerm =  pInfty_p^2 + pInfty_q^2 - (2 * scalarProd * pInfty_p * pInfty_q);
                coeff_q = 2*( (scalarProd * pInfty_p) - pInfty_q );
                coeff_p = 2*( (scalarProd * pInfty_q) - pInfty_p );
                coeff_pq = - 2 * scalarProd;
    
                e_1 = squareOperatorBasis(M+1+trI.numUsedEdges, 1, 1);
                e_p = squareOperatorBasis(M+1+trI.numUsedEdges, 1, ii+1);
                e_q = squareOperatorBasis(M+1+trI.numUsedEdges, 1, iN+1);
                e_pq = squareOperatorBasis(M+1+trI.numUsedEdges, ii+1, iN+1);
                e_p2 = squareOperatorBasis(M+1+trI.numUsedEdges, ii+1, ii+1);
                e_q2 = squareOperatorBasis(M+1+trI.numUsedEdges, iN+1, iN+1);
    
                e = (e_1 * constTerm) + (e_p * coeff_p) + (e_q * coeff_q) + (e_pq * coeff_pq) + e_p2 + e_q2;
                
                operatorList{iImg, ii, iNeighbor} = e;
    
                if(iNeighbor < trI.numMaxNeighbors)
                    cN = nng(ii, iNeighbor+1);
                    
                    [triangleIndex, edgeIndicesOrg] = trI.getIndex(ii, iN, cN);
    
                    if(triangleIndex > 0)
                        edgeIndices = zeros(1,3);
                        for ll = 1:3
                            edgeIndices(1,ll) = trI.getUtilisedEdgeIndex(edgeIndicesOrg(ll));
                        end
                        pInfty_r = norm(odReconstruction{iImg}(cN,:)) + EPSILON;
                        [G] = getQuadraticCoefficientsTriangleArea(directionVectorsList{iImg}(:,ii), directionVectorsList{iImg}(:,iN), ...
                                                                                        directionVectorsList{iImg}(:,cN));
                        [cg] = getCoefficientsTriangleAreaComplDepth(G, pInfty_p, pInfty_q, pInfty_r);
    
                        triVert = sort([ii iN cN]);
        
                        [eList] = getEquiarealOperatorBasis(triVert(1), triVert(2), triVert(3), edgeIndices, M, trI.numUsedEdges);

                        edgeIndex1 = edI.getIndex(ii, iN);
                        edgeIndex2 = edI.getIndex(ii, cN);
                        edgeIndex3 = edI.getIndex(iN, cN);
                        [eArea] = getHeronsOperatorBasis(edgeIndex1, edgeIndex2, edgeIndex3, edI.numberOfEdges);
                        areaOpBas = - eArea{1} + ( 2 .* eArea{2}) + ( 2 .* eArea{3}) - eArea{4} + ( 2 .* eArea{5}) -  eArea{6};
        
                        operatorListAreaFrmGeod{iImg,triangleIndex} = areaOpBas;
    
                        % Appending -ve area squared to operator basis's element (1,1)
                        cg(1) = cg(1);
    
                        eEqui = 0;
                        for iOp = 1:23
                            eEqui = eEqui + (cg(iOp)*eList{iOp});
                        end
                        
                        % If the same triangle is being visited multiple times,
                        % only the last pass will be stored, overwriting the
                        % previous operators. Thereby equiareal constraint will
                        % be applied only once per triangle (ToDo: verify if
                        % this approach is correct)
                        operatorListEquiareal{iImg,triangleIndex} = eEqui;
    
                        for kk = 1:numel(edgeIndices)
                            edgeIndex = edgeIndices(kk); 
                            if (edgeVisits(iImg, edgeIndex) == 0)
                                edgeVisits(iImg, edgeIndex) = 1;
                                if(kk == 1)
                                    eMomConPositive = squareOperatorBasis(M+1+trI.numUsedEdges, triVert(1)+1, triVert(2)+1);
                                elseif(kk == 2)
                                    eMomConPositive = squareOperatorBasis(M+1+trI.numUsedEdges, triVert(1)+1, triVert(3)+1);
                                elseif(kk == 3)
                                    eMomConPositive = squareOperatorBasis(M+1+trI.numUsedEdges, triVert(2)+1, triVert(3)+1);    
                                end
                                eMomConNegative = squareOperatorBasis(M+1+trI.numUsedEdges, 1, edgeIndex+1);
                                
                                momConsList{iImg,edgeIndex} = eMomConPositive - eMomConNegative;
                            end
                        end
                    end
                end
    
            end
        end
    end

    precisionVal = 1.35e-06;
    cvx_begin
    cvx_precision low
    if strcmp(solver, 'mosek')
        cvx_solver mosek
    elseif strcmp(solver, 'sedumi')
        cvx_solver sedumi
    else
        cvx_solver sdpt3
    end
    
        variable C(N, M+1+trI.numUsedEdges, M+1+trI.numUsedEdges);
        variable geod(1+edI.numberOfEdges,1+edI.numberOfEdges) nonnegative;

        isoCost = 0; equiarealCost = 0; traceCost = 0;

        for iImg = 1:N
            cC = squeeze(C(iImg, :, :));

            % % Equiareal cost [START]
            for ii = 1:M
                for iNeighbor = 1:size(nng, 2)
                    iN = nng(ii, iNeighbor);
                    edgeIndex = edI.getIndex(ii, iN);
                    isoCost = isoCost + abs(trace( operatorList{iImg, ii, iNeighbor} * cC ) - geod(1, edgeIndex));
                    if(iNeighbor < trI.numMaxNeighbors)
                        cN = nng(ii, iNeighbor+1);
                        [triangleIndex, ~] = trI.getIndex(ii, iN, cN);
                        if(triangleIndex > 0)
                            area = trace(operatorListAreaFrmGeod{iImg,triangleIndex} * geod);
                            equiarealCost = equiarealCost + abs(trace( operatorListEquiareal{triangleIndex} * cC ) - area);
                        end
                    end
                end
    
                cC(1,ii+1) <= norm(odReconstruction{iImg}(ii,:)) + EPSILON;
            end
            % % Equiareal cost [END]

            % % Pseudo-localizing constraints [START]
            for ii = (M+1+1):(M+1+trI.numUsedEdges)
            A = cC(1, ii);
                for jj = ii:(M+1+trI.numUsedEdges)
                    B = cC(jj, ii);
                    B <= A;
                end
            end
            % % Pseudo-localizing constraints [END]

            % % Moment constraints [START]
            for ii = 1:size(momConsList,2)
                if(numel(momConsList{iImg, ii}) > 0)
                    trace( momConsList{iImg, ii} * cC ) == 0;
                end
            end
            % % Moment constraints [END]            

            cC(1,1) == 1;
            cC == semidefinite(M+1+trI.numUsedEdges);

            traceCost = traceCost + trace(cC);
        end

        sum(geod(1,2:end)) == sumSqGeod;
        geod == semidefinite(1+edI.numberOfEdges);

        minimize( lambda*isoCost + lambda_equi*equiarealCost + (0.9*traceCost) + trace(geod));
    cvx_end


    rankList = [];
    for iImg = 1:N
        cC = squeeze(C(iImg, :, :));
        rankVal = effectiveRank(full(cC), precisionVal);
        gramMatrixList(iImg).G = full(cC);
        rankList = [rankList rankVal];
    
        if(rankVal~=1)
            [U, S, V] = svd(full(cC));
            D = U(:,1)*S(1,1)*V(:,1).';
            depthOfPoints = D(1,2:M+1);
        else
            depthOfPoints = cC(1,2:M+1);
        end
    
        for ii = 1:M
            depthOfPoints(ii) = norm(odReconstruction{iImg}(ii,:)) + EPSILON - depthOfPoints(ii);
        end
    
        reconstructionC = []; inextError = [];
        for ii = 1:M
            dirV = directionVectorsList{iImg}(:,ii).';
            P = dirV.*depthOfPoints(ii);
            reconstructionC = [reconstructionC; P];
    
            for iNeighbor = 1:size(nng, 2)
                    iN = nng(ii, iNeighbor);
                    dirVq = directionVectorsList{iImg}(:,iN).';
                    Q = dirVq.*depthOfPoints(iN);
                    edgeIndex = edI.getIndex(ii, iN);
                    inextError = [inextError abs(geod(edgeIndex) - norm(P - Q, 2)^2)];     
                end
        end
        fprintf('<strong>(Rank, trace, isoCost)</strong> of/from Gram matrix no. %d: = (%d, %d, %d)\n', iImg, rankVal, trace(full(cC)), isoCost);
        [regs] = absor(reconstructionC.', Data.Pgth(iImg).P, 'doScale', true, 'doTrans', false);
        reconsScaled = regs.s.*reconstructionC;
        reconstruction{iImg} = reconsScaled;
    end
    solutionStatus.optVal = cvx_optval;
    solutionStatus.cvx_status = cvx_status;
    solutionStatus.rank = rankList;        
    solutionStatus.gram_matrices = gramMatrixList;
    close all;

end

function [E] = squareOperatorBasis(M,i, j)
    E = zeros(M,M);
    E(i,j) = 1;
end

function [eList] = getEquiarealOperatorBasis(j, q, r, edges, M, numEdges)
    % % Edges follow the convention: [j-q, j-r, r-q]
    jq = M + edges(1); jr = M + edges(2); rq = M + edges(3);
    M = M + numEdges + 1;
    e{1} = squareOperatorBasis(M,1, 1);
    e{2} = squareOperatorBasis(M,jq+1, rq+1);
    e{3} = squareOperatorBasis(M,jr+1, rq+1);
    e{4} = squareOperatorBasis(M,jq+1, r+1);
    e{5} = squareOperatorBasis(M,jq+1, jr+1);
    e{6} = squareOperatorBasis(M,rq+1, rq+1);
    e{7} = squareOperatorBasis(M,jq+1, jq+1);
    e{8} = squareOperatorBasis(M,jr+1, jr+1);
    e{9} = squareOperatorBasis(M,r+1, j+1);
    e{10} = squareOperatorBasis(M,q+1, j+1);
    e{11} = squareOperatorBasis(M,jq+1, j+1);
    e{12} = squareOperatorBasis(M,q+1, jq+1);
    e{13} = squareOperatorBasis(M,jr+1, j+1);
    e{14} = squareOperatorBasis(M,r+1, jr+1);
    e{15} = squareOperatorBasis(M,rq + 1, q+1);
    e{16} = squareOperatorBasis(M,r+1, rq+1);
    e{17} = squareOperatorBasis(M,j+1, j+1);
    e{18} = squareOperatorBasis(M,r+1, q+1);
    e{19} = squareOperatorBasis(M,r+1, r+1);
    e{20} = squareOperatorBasis(M,q+1, q+1);
    e{21} = squareOperatorBasis(M,1, j+1);
    e{22} = squareOperatorBasis(M,1, q+1);
    e{23} = squareOperatorBasis(M,1, r+1);

    eList = e;
end

function [eList] = getHeronsOperatorBasis(e1, e2, e3, numEdges)
    % Equation: -a^4 + 2*a^2*b^2 + 2*a^2*c^2 - b^4 + 2*b^2*c^2 - c^4
    % Equation after change of variable: -a^2 + 2*a*b + 2*a*c - b^2 + 2*b*c - c^2

    c{1} = squareOperatorBasis(1+numEdges, e1, e1);
    c{2} = squareOperatorBasis(1+numEdges, e1, e2);
    c{3} = squareOperatorBasis(1+numEdges, e1, e3);
    c{4} = squareOperatorBasis(1+numEdges, e2, e2);
    c{5} = squareOperatorBasis(1+numEdges, e2, e3);
    c{6} = squareOperatorBasis(1+numEdges, e3, e3);

    eList = c;
end

