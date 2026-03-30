function [reconstruction, solutionStatus] = SfT_IsometricZeroth(template, pts2d, K, nng, pNormType, lambda, precisionVal, solver)


    if ~exist('pNormType', 'var')
        pNormType = 'SqMOD';
    else
        assert(strcmp(pNormType, 'SqMOD') || strcmp(pNormType, 'MOD'), 'The only options for pNormType are <MOD> and <SqMOD> [default]');
    end    

    if ~exist('lambda', 'var')
        lambda = 1;
    end    

    if ~exist('precisionVal', 'var')
        precisionVal = 1.35e-06;
    else
        assert( (precisionVal < 1) && (precisionVal > 0), 'Precision value must be in the range [0, 1]' );
    end

    if ~exist('solver', 'var')
        solver = 'sdpt3';
    end    

    M = size(template,1);
    Kinv = inv(K); 
    directionVectors = [];

    for iPts = 1:M
        dirVec_j = Kinv*[pts2d(iPts,1); pts2d(iPts,2); 1];
        dirVec_j = dirVec_j./norm(dirVec_j);
        directionVectors = [directionVectors dirVec_j];
    end
    
    fprintf('\n\n ------- Isometric SfT [<strong>forward-depths</strong>] -------\n\n');

    cvx_begin
    cvx_precision high

    if strcmp(solver, 'mosek')
        cvx_solver mosek
    elseif strcmp(solver, 'sedumi')
        cvx_solver sedumi
    else
        cvx_solver sdpt3
    end    
    
        variable C(1+M,1+M);

        isoCost = 0; mdhCost = 0;

        for ii = 1:M
            p = template(ii,1:2);

            for iNeighbor = 1:size(nng, 2)
                iN = nng(iPts, iNeighbor);
                q = template(iN,1:2);
                scalarProd = dot(directionVectors(:,ii), directionVectors(:,iN));
                isoCost = isoCost + abs(C(ii+1,ii+1) + C(iN+1,iN+1) - (2*C(ii+1,iN+1)*scalarProd) - norm(p - q, 2)^2);
            end

            C(ii+1,ii+1) >= 0;

            if(strcmp(pNormType, 'MOD'))
                mdhCost = mdhCost - (C(1,ii+1));
            else
                mdhCost = mdhCost - (C(1,ii+1));
            end
        end

        C(1,1) == 1;
        C == semidefinite(M+1);


        minimize( (10*isoCost) + (0.1*trace(C))  + lambda*mdhCost);
    
    cvx_end

    rankVal = effectiveRank(full(C), precisionVal);

    solutionStatus.optVal = cvx_optval;
    solutionStatus.cvx_status = cvx_status;
    solutionStatus.rank = rankVal;
    solutionStatus.PreProjectionGramMatrix = full(C);

    if(rankVal~=1)
        [U, S, V] = svd(full(C));
        D = U(:,1)*S(1,1)*V(:,1).';
        depthOfPoints = D(1,2:M+1);
    else
        depthOfPoints = C(1,2:M+1);
    end

    reconstruction = []; inextError = [];
    for ii = 1:M
        p = template(ii,1:2);
        dirV = directionVectors(:,ii).';
        P = dirV.*depthOfPoints(ii);
        reconstruction = [reconstruction; P];

        for iNeighbor = 1:size(nng, 2)
                iN = nng(iPts, iNeighbor);
                q = template(iN,1:2);
                dirVq = directionVectors(:,iN).';
                Q = dirVq.*depthOfPoints(iN);
                inextError = [inextError abs(norm(p-q,2) - norm(P - Q, 2))];     
            end
    end

end