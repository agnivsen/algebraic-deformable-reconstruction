function [reconstruction, solutionStatus] = algSfT_mEqui(template, pts2d, K, nng, eps, solver, lambda, mode, precisionVal, priorRecons)

    if ~exist('eps', 'var')
        EPSILON = 0.1;
    else
        EPSILON = eps;
    end
    
    if ~exist('solver', 'var')
        solver = 'sdpt3';
    end

    % lambda{1} is weight for isometry, 
    % lambda{2} is weight for trace, 
    % lambda{3} is weight for optional costs (MDH, conformal, equiareal, etc.)
    if ~exist('lambda', 'var')
        lambda{1} = 10;
        lambda{2} = 0.001;
        lambda{3} = lambda{1}/10;
    end    

    if ~exist('mode', 'var')
        mode = 'approximate';
    else
        assert(strcmp(mode, 'approximate') || strcmp(mode, 'exact'), 'The only options for MODE are <approximate> [default] and <exact>');
    end        

    M = size(template,1);
    Kinv = inv(K); 
    directionVectors = [];

    if ~exist('priorRecons', 'var')
        pNormType = 'SqMOD';
        [reconstructionPriorFD, solStatFD] = SfT_IsometricZeroth(template, pts2d, K, nng, pNormType, lambda{4}, precisionVal, solver);
    else
        reconstructionPriorFD = priorRecons;
        solStatFD = false;
    end

    if ~exist('precisionVal', 'var')
        precisionVal = 1.35e-06;
    end      

    for iPts = 1:M
        dirVec_j = Kinv*[pts2d(iPts,1); pts2d(iPts,2); 1];
        dirVec_j = dirVec_j./norm(dirVec_j);
        directionVectors = [directionVectors dirVec_j];
    end

    trI = TriangleIndexer(nng);

    assert(trI.numberOfTriangles > 1, 'There has to be at least 1 valid triangle in the NNG to impose equiareal constraint!');

    fprintf(['\n(Equiareal SfT has <strong>%d triangles</strong> with <strong>%d total</strong> edges,... \n<strong>%d utilised</strong> edges and... \n<strong>moment matrix</strong> of size <strong>[%d x %d]</strong>' ...
        ' in this data).\n\n'], trI.numberOfTriangles, trI.numberOfEdges, trI.numUsedEdges, M+1+trI.numUsedEdges, M+1+trI.numUsedEdges);
    
    edgeVisits = zeros(1, trI.numUsedEdges);

    for ii = 1:M
        p = template(ii,:);
        pInfty_p = norm(reconstructionPriorFD(ii,:)) + EPSILON;

        for iNeighbor = 1:size(nng, 2)
            iN = nng(ii, iNeighbor);
            q = template(iN,:);
            pInfty_q = norm(reconstructionPriorFD(iN,:)) + EPSILON;
            scalarProd = dot(directionVectors(:,ii), directionVectors(:,iN));

            constTerm =  pInfty_p^2 + pInfty_q^2 - (2 * scalarProd * pInfty_p * pInfty_q) - norm(p - q, 2)^2 ;
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
            
            operatorList{ii, iNeighbor} = e;

            if(iNeighbor < trI.numMaxNeighbors)
                cN = nng(ii, iNeighbor+1);
                
                r = template(cN,:);

                [triangleIndex, edgeIndicesOrg] = trI.getIndex(ii, iN, cN);

                if(triangleIndex > 0)
                    edgeIndices = zeros(1,3);
                    for ll = 1:3
                        edgeIndices(1,ll) = trI.getUtilisedEdgeIndex(edgeIndicesOrg(ll));
                    end
                    pInfty_r = norm(reconstructionPriorFD(cN,:)) + EPSILON;
                    [G] = getQuadraticCoefficientsTriangleArea(directionVectors(:,ii), directionVectors(:,iN), directionVectors(:,cN));
                    [cg] = getCoefficientsTriangleAreaComplDepth(G, pInfty_p, pInfty_q, pInfty_r);

                    triVert = sort([ii iN cN]);
    
                    [eList] = getEquiarealOperatorBasis(triVert(1), triVert(2), triVert(3), edgeIndices, M, trI.numUsedEdges);

                    [area] = heronsFormula(norm(p - q, 2), norm(p - r, 2), norm(r - q, 2));

                    cg(1) = cg(1) - area^2;

                    eEqui = 0;
                    for iOp = 1:23
                        eEqui = eEqui + (cg(iOp)*eList{iOp});
                    end
                    
                    % If the same triangle is being visited multiple times,
                    % only the last pass will be stored, overwriting the
                    % previous operators. Thereby equiareal constraint will
                    % be applied only once per triangle 
                    operatorListEquiareal{triangleIndex} = eEqui;

                    for kk = 1:numel(edgeIndices)
                        edgeIndex = edgeIndices(kk); 
                        if (edgeVisits(edgeIndex) == 0)
                            edgeVisits(edgeIndex) = 1;
                            if(kk == 1)
                                eMomConPositive = squareOperatorBasis(M+1+trI.numUsedEdges, triVert(1)+1, triVert(2)+1);
                            elseif(kk == 2)
                                eMomConPositive = squareOperatorBasis(M+1+trI.numUsedEdges, triVert(1)+1, triVert(3)+1);
                            elseif(kk == 3)
                                eMomConPositive = squareOperatorBasis(M+1+trI.numUsedEdges, triVert(2)+1, triVert(3)+1);    
                            end
                            eMomConNegative = squareOperatorBasis(M+1+trI.numUsedEdges, 1, edgeIndex+M+1);
                            
                            momConsList{edgeIndex} = eMomConPositive - eMomConNegative;
                        end
                    end
                end
            end

        end
    end

    fprintf('\n\n ------- Equiareal SfT [<strong>opposite-depths</strong>] -------\n\n');

    cvx_begin
    cvx_precision low
    if strcmp(solver, 'mosek')
        cvx_solver mosek
    elseif strcmp(solver, 'sedumi')
        cvx_solver sedumi
    else
        cvx_solver sdpt3
    end
    
        variable C(M+1+trI.numUsedEdges, M+1+trI.numUsedEdges);

        isoCost = 0; equiarealCost = 0;

        for ii = 1:M
            for iNeighbor = 1:size(nng, 2)
                iN = nng(ii, iNeighbor);
                isoCost = isoCost + abs(trace( operatorList{ii, iNeighbor} * C ));
                if(iNeighbor < trI.numMaxNeighbors)
                    cN = nng(ii, iNeighbor+1);
                    [triangleIndex, ~] = trI.getIndex(ii, iN, cN);
                    if(triangleIndex > 0)
                        equiarealCost = equiarealCost + abs(trace( operatorListEquiareal{triangleIndex} * C ));
                    end
                end
            end

            C(1,ii+1) <= norm(reconstructionPriorFD(ii,:)) + EPSILON;
        end


        % % Moment constraints [START]
        for ii = 1:size(momConsList,2)
            if(numel(momConsList{ii}) > 0)
                trace( momConsList{ii} * C ) == 0;
            end
        end
        % % Moment constraints [END]

        C(1,1) == 1;
        C == semidefinite(M+1+trI.numUsedEdges);

        minimize( (lambda{1}*isoCost) + (lambda{2} * trace(C)) + (lambda{3}*equiarealCost) );
    cvx_end

    rankVal = effectiveRank(full(C), precisionVal);

    solutionStatus.optVal = cvx_optval;
    solutionStatus.cvx_status = cvx_status;
    solutionStatus.rank = rankVal;    
    solutionStatus.PreProjectionGramMatrix = full(C);

    solutionStatus.FDrecons = reconstructionPriorFD;
    solutionStatus.FDsolStat = solStatFD;    

    if(rankVal~=1)
        [U, S, V] = svd(full(C));
        D = U(:,1)*S(1,1)*V(:,1).';
        depthOfPoints = D(1,2:M+1);
    else
        depthOfPoints = C(1,2:M+1);
    end

    for ii = 1:M
        depthOfPoints(ii) = norm(reconstructionPriorFD(ii,:)) + EPSILON - depthOfPoints(ii);
    end

    reconstruction = []; inextError = [];
    for ii = 1:M
        p = template(ii,:);
        dirV = directionVectors(:,ii).';
        P = dirV.*depthOfPoints(ii);
        reconstruction = [reconstruction; P];

        for iNeighbor = 1:size(nng, 2)
                iN = nng(iPts, iNeighbor);
                q = template(iN,:);
                dirVq = directionVectors(:,iN).';
                Q = dirVq.*depthOfPoints(iN);
                inextError = [inextError abs(norm(p-q,2) - norm(P - Q, 2))];     
            end
    end

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

