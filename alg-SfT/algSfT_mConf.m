function [reconstruction, solutionStatus] = algSfT_mConf(template, pts2d, K, nng, eps, solver, lambda, mode, precisionVal, priorRecons)

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
        lambda{3} = lambda{1};
    end    

    if ~exist('mode', 'var')
        mode = 'approximate';
    else
        assert(strcmp(mode, 'approximate') || strcmp(mode, 'exact'), 'The only options for MODE are <approximate> [default] and <exact>');
    end    

    if ~exist('precisionVal', 'var')
        precisionVal = 1.35e-06;
    else
        assert( (precisionVal < 1) && (precisionVal > 0), 'Precision value must be in the range [0, 1]' );
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

    for iPts = 1:M
        dirVec_j = Kinv*[pts2d(iPts,1); pts2d(iPts,2); 1];
        dirVec_j = dirVec_j./norm(dirVec_j);
        directionVectors = [directionVectors dirVec_j];
    end
    

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

            e_1 = squareOperatorBasis(M+1, 1, 1);
            e_p = squareOperatorBasis(M+1, 1, ii+1);
            e_q = squareOperatorBasis(M+1, 1, iN+1);
            e_pq = squareOperatorBasis(M+1, ii+1, iN+1);
            e_p2 = squareOperatorBasis(M+1, ii+1, ii+1);
            e_q2 = squareOperatorBasis(M+1, iN+1, iN+1);

            e = (e_1 * constTerm) + (e_p * coeff_p) + (e_q * coeff_q) + (e_pq * coeff_pq) + e_p2 + e_q2;
            
            operatorList{ii, iNeighbor} = e;

            if(iNeighbor < size(nng, 2))
                cN = nng(ii, iNeighbor+1);
                r = template(cN,:);
                pInfty_r = norm(reconstructionPriorFD(cN,:)) + EPSILON;
                scalarProd_r = dot(directionVectors(:,ii), directionVectors(:,cN));
                [confCoeff] = getConfCoeff(pInfty_p, pInfty_q, pInfty_r, scalarProd_r, scalarProd, norm(p - q, 2)^2, norm(p - r, 2)^2);

                e_r = squareOperatorBasis(M+1, 1, cN+1);
                e_r2 = squareOperatorBasis(M+1, cN+1, cN+1);
                e_pr = squareOperatorBasis(M+1, ii+1, cN+1);

                eConf = ( confCoeff(1) * e_p2) + ( confCoeff(2) * e_q2) + ( confCoeff(3) * e_r2) + ( confCoeff(4) * e_p) + ...
                                        ( confCoeff(5) * e_q) + ( confCoeff(6) * e_r) + ( confCoeff(7) * e_pq) + ( confCoeff(8) * e_pr) + ...
                                            ( confCoeff(9) * e_1);
                operatorListConformal{ii, iNeighbor} = eConf;
            end

        end
    end

    fprintf('\n\n ------- Conformal SfT [<strong>opposite-depths</strong>] -------\n\n');

    cvx_begin
    cvx_precision low
    if strcmp(solver, 'mosek')
        cvx_solver mosek
    elseif strcmp(solver, 'sedumi')
        cvx_solver sedumi
    else
        cvx_solver sdpt3
    end
    
        variable C(M+1,M+1);

        isoCost = 0; confCost = 0;

        for ii = 1:M
            for iNeighbor = 1:size(nng, 2)
                iN = nng(ii, iNeighbor);
                isoCost = isoCost + abs(trace( operatorList{ii, iNeighbor} * C ));
                if(iNeighbor < size(nng, 2))
                    confCost = confCost + abs(trace( operatorListConformal{ii, iNeighbor} * C ));
                end
            end

            C(1,ii+1) <= norm(reconstructionPriorFD(ii,:)) + EPSILON;
        end


        C(1,1) == 1;
        C == semidefinite(M+1);

        minimize( (lambda{1}*isoCost) + (lambda{2} * trace(C)) + (lambda{3}*confCost) );
    
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

function [conformalCoeff] = getConfCoeff(infty_j, infty_q, infty_r, alpha_jr, alpha_jq, d_jq, d_jr)
    constTerm = -2 * infty_j * infty_q * alpha_jq * d_jr + 2 * infty_j * infty_r * alpha_jr * d_jq - infty_j ^ 2 * d_jq + infty_j ^ 2 * d_jr + infty_q ^ 2 * d_jr - infty_r ^ 2 * d_jq;
    
    coeff_j_sqr = d_jr - d_jq;
    
    coeff_q_sqr = d_jr;
    
    coeff_r_sqr = -d_jq;
    
    coeff_j = 2 * infty_q * alpha_jq * d_jr - 2 * infty_r * alpha_jr * d_jq + 2 * infty_j * d_jq - 2 * infty_j * d_jr;
    
    coeff_q =2 * infty_j * alpha_jq * d_jr - 2 * infty_q * d_jr;
    
    coeff_r = -2 * infty_j * alpha_jr * d_jq + 2 * infty_r * d_jq;
    
    coeff_jq = -2*alpha_jq*d_jr;
    
    coeff_jr = 2*alpha_jr*d_jq;

    conformalCoeff = [coeff_j_sqr coeff_q_sqr coeff_r_sqr coeff_j coeff_q coeff_r coeff_jq coeff_jr constTerm];
end
