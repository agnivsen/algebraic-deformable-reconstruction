function [reconstruction, solutionStatus] = algSfT_mIso(template, pts2d, K, nng, eps, solver, lambda, mode, precisionVal, pNormType, variant, priorRecons)

    if ~exist('eps', 'var')
        EPSILON = 0.1;
    else
        EPSILON = eps;
    end
    
    if ~exist('solver', 'var')
        solver = 'sdpt3';
    end

    if ~exist('mode', 'var')
        mode = 'approximate';
    else
        assert(strcmp(mode, 'approximate') || strcmp(mode, 'exact'), 'The only options for MODE are <approximate> [default] and <exact>');
    end

    if ~exist('pNormType', 'var')
        pNormType = 'SqMOD';
    else
        assert(strcmp(pNormType, 'SqMOD') || strcmp(pNormType, 'MOD'), 'The only options for pNormType are <MOD> and <SqMOD> [default]');
    end    

    if ~exist('variant', 'var')
        variant = 'Deformable';
    else
        assert(strcmp(variant, 'Deformable') || strcmp(variant, 'Rigid'), 'The only options for variant are <Rigid> and <Deformable> [default]');
    end    

    if ~exist('precisionVal', 'var')
        precisionVal = 1.35e-06;
    else
        assert( (precisionVal < 1) && (precisionVal > 0), 'Precision value must be in the range [0, 1]' );
    end

    %

    % lambda{1} is weight for isometry, 
    % lambda{2} is weight for trace, 
    % lambda{3} is weight for optional costs (MDH, conformal, equiareal, etc.)
    if ~exist('lambda', 'var')
        lambda{1} = 10;
        lambda{2} = 0.001;
        lambda{3} = 1;
        lambda{4} = 1;
    end    

    M = size(template,1);
    Kinv = inv(K);  
    directionVectors = [];


    if ~exist('priorRecons', 'var')
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

        end
    end

    fprintf('\n\n ------- Isometric SfT [<strong>opposite-depths</strong>] -------\n\n');

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

        isoCost = 0;
        if(strcmp(pNormType, 'MOD'))
            mdhCost = 0;
        end

        for ii = 1:M
            for iNeighbor = 1:size(nng, 2)
                if(strcmp(variant, 'Rigid'))
                    trace( operatorList{ii, iNeighbor} * C) ==0;
                else
                    isoCost = isoCost + abs(trace( operatorList{ii, iNeighbor} * C ));
                end                
            end

            if(strcmp(pNormType, 'MOD'))
                mdhCost = mdhCost + abs(C(1, ii+1));
            end
        end


        C(1,1) == 1;
        C == semidefinite(M+1);

        if(strcmp(pNormType, 'MOD'))
            minimize( (lambda{1}*isoCost) + lambda{2}*trace(C) + lambda{3}*mdhCost);
        else
            minimize( (lambda{1}*isoCost) + lambda{2}*trace(C) );
        end
    
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
