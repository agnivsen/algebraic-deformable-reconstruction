function [reconstruction, solutionStatus] = algNRSfM_mConf(Data, K, nng, odReconstruction, sumSqGeod, eps, solver, lambda)

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

    for iImg = 1:N
        directionVectors = [];
        for iPts = 1:M
            dirVec_j = Kinv*[Data.p(iImg).p(1:2, iPts); 1];
            dirVec_j = dirVec_j./norm(dirVec_j);
            directionVectors = [directionVectors dirVec_j];
        end
        directionVectorsList{iImg} = directionVectors;
    end
    
    lambda_conf = lambda/100;

    % Gram matrix from vector: [1, delta_1, delta_2, ..., delta_m, delta_1^2, delta_2^2, ..., delta_m^2, delta_1 delta_2, .... , edge_1, edle_2, ..., edge_k]

    gramMatrixDim = 1 + (2*M) + (2*edI.numberOfEdges);

    for iImg = 1:N
        for ii = 1:M
            pInfty_p = norm(odReconstruction{iImg}(ii,:)) + EPSILON;

            momConDeltaPositive = squareOperatorBasis(gramMatrixDim, 1+ii, 1+ii);
            momConDeltaNegative = squareOperatorBasis(gramMatrixDim, 1, 1+M+ii);
            momConDelta{ii} = momConDeltaPositive - momConDeltaNegative;
    
            for iNeighbor = 1:size(nng, 2)
                iN = nng(ii, iNeighbor);
                pInfty_q = norm(odReconstruction{iImg}(iN,:)) + EPSILON;
                scalarProd = dot(directionVectorsList{iImg}(:,ii), directionVectorsList{iImg}(:,iN));

                constTerm =  pInfty_p^2 + pInfty_q^2 - (2 * scalarProd * pInfty_p * pInfty_q) ;
                coeff_q = 2*( (scalarProd * pInfty_p) - pInfty_q );
                coeff_p = 2*( (scalarProd * pInfty_q) - pInfty_p );
                coeff_pq = - 2 * scalarProd;
    
                e_1 = squareOperatorBasis(gramMatrixDim, 1, 1);
                e_p = squareOperatorBasis(gramMatrixDim, 1, ii+1);
                e_q = squareOperatorBasis(gramMatrixDim, 1, iN+1);
                e_pq = squareOperatorBasis(gramMatrixDim, ii+1, iN+1);
                e_p2 = squareOperatorBasis(gramMatrixDim, ii+1, ii+1);
                e_q2 = squareOperatorBasis(gramMatrixDim, iN+1, iN+1);
    
                e = (e_1 * constTerm) + (e_p * coeff_p) + (e_q * coeff_q) + (e_pq * coeff_pq) + e_p2 + e_q2;
                
                operatorList{ii, iNeighbor} = e;

                edgeIndex_jq = edI.getIndex(ii, iN);

                momConEdgePositive = squareOperatorBasis(gramMatrixDim, 1+ii, 1+iN);
                momConEdgeNegative = squareOperatorBasis(gramMatrixDim, 1, 1+(2*M) + edgeIndex_jq);
                momConEdge{edgeIndex_jq} = momConEdgePositive - momConEdgeNegative;
    
                if(iNeighbor < size(nng, 2))
                    cN = nng(ii, iNeighbor+1);
                    pInfty_r = norm(odReconstruction{iImg}(cN,:)) + EPSILON;
                    scalarProd_r = dot(directionVectors(:,ii), directionVectors(:,cN));
                    %  coeff_j_sqr_d_jr,  coeff_j_sqr_d_jq,  coeff_q_sqr_d_jr,  coeff_r_sqr_d_jq,  coeff_j_d_jr,  coeff_j_d_jq, ...
                    %     coeff_q_d_jr,  coeff_r_d_jq,  coeff_jq_d_jr,  coeff_jr_d_jq,  coeff_d_jr,  coeff_d_jq
                    [confCoeff] = getConfCoeff(pInfty_p, pInfty_q, pInfty_r, scalarProd_r, scalarProd);

                    edgeIndex_jr = edI.getIndex(ii, cN);

                    Ij = 1+ii; Iq = 1+iN; Ir = 1+cN;    
                    Ij_sqr = 1 + M + ii; Iq_sqr = 1 + M + iN; Ir_sqr = 1 + M + cN;
                    Ijq = 1 + (2*M) + edgeIndex_jq; Ijr = 1 + (2*M) + edgeIndex_jr;
                    Id_jq = 1 + (2*M) + edI.numberOfEdges + edgeIndex_jq; Id_jr = 1 + (2*M) + edI.numberOfEdges + edgeIndex_jr;
    
                    eC{1} = squareOperatorBasis(gramMatrixDim, Ij_sqr, Id_jr);
                    eC{2} = squareOperatorBasis(gramMatrixDim, Ij_sqr, Id_jq);
                    eC{3} = squareOperatorBasis(gramMatrixDim, Iq_sqr, Id_jr);
                    eC{4} = squareOperatorBasis(gramMatrixDim, Ir_sqr, Id_jq);
                    eC{5} = squareOperatorBasis(gramMatrixDim, Ij, Id_jr);
                    eC{6} = squareOperatorBasis(gramMatrixDim, Ij, Id_jq);
                    eC{7} = squareOperatorBasis(gramMatrixDim, Iq, Id_jr);
                    eC{8} = squareOperatorBasis(gramMatrixDim, Ir, Id_jq);
                    eC{9} = squareOperatorBasis(gramMatrixDim, Ijq, Id_jr);
                    eC{10} = squareOperatorBasis(gramMatrixDim, Ijr, Id_jq);
                    eC{11} = squareOperatorBasis(gramMatrixDim, 1, Id_jr);
                    eC{12} = squareOperatorBasis(gramMatrixDim, 1, Id_jq);
    
                    eConf = 0;
                    for iCoeff = 1:12
                        eConf = eConf + (confCoeff(iCoeff).*eC{iCoeff});
                    end
                    operatorListConformal{ii, iNeighbor} = eConf;
                end
    
            end
        end
    end

    precisionVal = 1.35e-06;
    fprintf('---: Conformal NRSfM: over to CVX\n\n');
    cvx_begin
    cvx_precision low
    if strcmp(solver, 'mosek')
        cvx_solver mosek
    elseif strcmp(solver, 'sedumi')
        cvx_solver sedumi
    else
        cvx_solver sdpt3
    end
    
        variable C(N, gramMatrixDim, gramMatrixDim);

        isoCost = 0; confCost = 0; traceCost = 0;

        geod_1 = squeeze(C(1, 1, 1 + (2*M) + (edI.numberOfEdges)+1:end));
        for iImg = 1:N
            fprintf('==: (NRSfM conformal) : modelling image #<strong>%d</strong> of %d\n', iImg, N);
            for ii = 1:M
                cC = squeeze(C(iImg, :, :));
                for iNeighbor = 1:size(nng, 2)
                    iN = nng(ii, iNeighbor);
                    edgeIndex_jq = edI.getIndex(ii, iN);
                    isoCost = isoCost + abs(trace( operatorList{ii, iNeighbor} * cC ) - cC(1,1 + (2*M) + (edI.numberOfEdges)+edgeIndex_jq));                
                    trace(momConEdge{edgeIndex_jq}*cC) == 0;
                    if(iNeighbor < size(nng, 2))
                        confCost = confCost + abs(trace( operatorListConformal{ii, iNeighbor} * cC ));
                    end
                end
    
                trace(momConDelta{ii}*cC) == 0;
    
                cC(1,ii+1) <= norm(odReconstruction{iImg}(ii,:)) + EPSILON;
                cC(1,1) == 1;
                cC == semidefinite(gramMatrixDim);
                sum(cC(1 + (2*M) + (edI.numberOfEdges)+1:end, 1 + (2*M) + (edI.numberOfEdges)+1:end)) == sumSqGeod;
                traceCost = traceCost + trace(cC);
            end

            % Tensor constraint (for geodesics)
            if(iImg == 1)
                for iImg2 = 2:N
                    geod_n = squeeze(C(iImg2, 1, 1 + (2*M) + (edI.numberOfEdges)+1:end));
                    geod_1 == geod_n;
                end
            end
        end

        fprintf('==: (NRSfM conformal) : done modelling\n');

        minimize( lambda*isoCost + lambda_conf*confCost + (0.9*traceCost) );
    
    cvx_end

    rankList = [];
    for iImg = 1:N
        cC = squeeze(C(iImg, :, :));
        rankVal = effectiveRank(full(cC), precisionVal);
        rankList = [rankList rankVal];
        gramMatrixList(iImg).G = full(cC);
        geod = cC(1 + (2*M) + (edI.numberOfEdges)+1:end, 1 + (2*M) + (edI.numberOfEdges)+1:end);
    
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

function [conformalCoeff] = getConfCoeff(infty_j, infty_q, infty_r, alpha_jr, alpha_jq)

    coeff_j_sqr_d_jr = 1;
    coeff_j_sqr_d_jq = -1;
    coeff_q_sqr_d_jr = 1;
    coeff_r_sqr_d_jq = -1;
    coeff_j_d_jr = ((2 * infty_q * alpha_jq) - ( 2 * infty_j));
    coeff_j_d_jq = - ((2 * infty_r * alpha_jr) + (2 * infty_j));
    coeff_q_d_jr = (2 * infty_j * alpha_jq ) - (2 * infty_q );
    coeff_r_d_jq = -(2 * infty_j * alpha_jr) + (2 * infty_r);
    coeff_jq_d_jr = -2*alpha_jq;
    coeff_jr_d_jq = 2*alpha_jr;
    coeff_d_jr = ((-2 * infty_j * infty_q * alpha_jq)  + infty_j^2 + infty_q^2 );
    coeff_d_jq = - (infty_r^2 + (2 * infty_j * infty_r * alpha_jr ) - infty_j^2);

    conformalCoeff = [coeff_j_sqr_d_jr coeff_j_sqr_d_jq coeff_q_sqr_d_jr coeff_r_sqr_d_jq coeff_j_d_jr coeff_j_d_jq ...
        coeff_q_d_jr coeff_r_d_jq coeff_jq_d_jr coeff_jr_d_jq coeff_d_jr coeff_d_jq] ;
end
