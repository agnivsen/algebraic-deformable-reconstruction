function [reconstruction, solutionStatus] = algNRSfM_mIso(Data, K, nng, odReconstruction, sumSqGeod, eps, solver, lambda, pNormType)

    if ~exist('eps', 'var')
        EPSILON = 0.1;
    else
        EPSILON = eps;
    end
    
    if ~exist('solver', 'var')
        solver = 'sdpt3';
    end

    if ~exist('pNormType', 'var')
        normType = 'SqMOD';
    else
        assert(strcmp(pNormType, 'SqMOD') || strcmp(pNormType, 'MOD'), 'The only options for pNormType are <MOD> and <SqMOD>');
        normType = pNormType;
    end    

    if ~exist('lambda', 'var')
        lambda = 10;
    end

    N = size(Data.p,2);
    M = size(Data.p(1).p,2);
    Kinv = inv(K);  
    edI = EdgeIndexer(nng);

    for iImg = 1:N
        directionVectors = [];
        for iPts = 1:M
            dirVec_j = Kinv*[Data.p(iImg).p(1:2, iPts); 1];
            dirVec_j = dirVec_j./norm(dirVec_j);
            directionVectors = [directionVectors dirVec_j];
        end
        directionVectorsList{iImg} = directionVectors;
    end
    

    for iImg = 1:N
        for ii = 1:M
            pInfty_p = norm(odReconstruction{iImg}(ii,:)) + EPSILON;
    
            for iNeighbor = 1:size(nng, 2)
                iN = nng(ii, iNeighbor);
                pInfty_q = norm(odReconstruction{iImg}(iN,:)) + EPSILON;
                scalarProd = dot(directionVectorsList{iImg}(:,ii), directionVectorsList{iImg}(:,iN));
    
                constTerm =  pInfty_p^2 + pInfty_q^2 - (2 * scalarProd * pInfty_p * pInfty_q)  ;
                coeff_q = 2*( (scalarProd * pInfty_p) - pInfty_q );
                coeff_p = 2*( (scalarProd * pInfty_q) - pInfty_p );
                coeff_pq = - 2 * scalarProd;
    
                e_1 = squareOperatorBasis(M+1, 1, 1);
                e_p = squareOperatorBasis(M+1, 1, ii+1);
                e_q = squareOperatorBasis(M+1, 1, iN+1);
                e_pq = squareOperatorBasis(M+1, ii+1, iN+1);
                e_p2 = squareOperatorBasis(M+1, ii+1, ii+1);
                e_q2 = squareOperatorBasis(M+1, iN+1, iN+1);
    
                e = (e_1 * constTerm) + (e_p * coeff_p) + (e_q * coeff_q) + (e_pq * coeff_pq) + e_p2 + e_q2;
                
                operatorList{iImg, ii, iNeighbor} = e;
    
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
    
        variable C(N, M+1,M+1);
        variable geod(edI.numberOfEdges) nonnegative;

        isoCost = 0; traceCost = 0;
        if(strcmp(pNormType, 'MOD'))
            mdhCost = 0;
        end

        for iImg = 1:N
            cC = squeeze(C(iImg, :, :));
            for ii = 1:M
                for iNeighbor = 1:size(nng, 2)
                    iN = nng(ii, iNeighbor);
                    edgeIndex = edI.getIndex(ii, iN);
                    isoCost = isoCost + abs(trace( operatorList{iImg, ii, iNeighbor} * cC ) - geod(edgeIndex));
                end
    
                if(strcmp(pNormType, 'MOD'))
                    mdhCost = mdhCost + abs(cC(1, ii+1));
                end
    
                cC(1,ii+1) <= norm(odReconstruction{iImg}(ii,:)) + EPSILON;
            end
            cC(1,1) == 1;
            cC(:, :) == semidefinite(M+1);
            traceCost = traceCost + trace(cC);
        end

        sum(geod) == sumSqGeod;

        if(strcmp(pNormType, 'MOD'))
            minimize( (lambda*isoCost) + traceCost + mdhCost);
        else
            minimize( (lambda*isoCost) + 0.1*traceCost );
        end
    
    cvx_end

    rankList = []; 
    for iImg = 1:N
        cC = squeeze(C(iImg, :, :));
        rankVal = effectiveRank(full(cC), precisionVal);
        rankList = [rankList rankVal];
        gramMatrixList(iImg).G = full(cC);
    
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
        [regs] = absor(reconstructionC.', Data.Pgth(iImg).P, 'doScale', true, 'doTrans', false);
        reconsScaled = regs.s.*reconstructionC;
        reconstruction{iImg} = reconsScaled;
    end

    solutionStatus.optVal = cvx_optval;
    solutionStatus.cvx_status = cvx_status;
    solutionStatus.rank = rankList;    
    solutionStatus.gram_matrices = gramMatrixList;

end

function [E] = squareOperatorBasis(M,i, j)
    E = zeros(M,M);
    E(i,j) = 1;
end
