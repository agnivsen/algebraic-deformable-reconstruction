% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function pp = paperComputeBe(pp,sh,q_FP,Q_FP)

altot = sort(pp.al(:)');
% bending angle initialisation
% distance computation
[d,alnth] = newPaperComputeD(pp,sh);

% bending angle computation
ns = 50;

% TPS
E = GTPSW_TPS_compute_FDPMat(q_FP');

PE = Q_FP*E';
thn = 1e-7;
dtest = 1e-3;

for i=1:length(pp.be)
    % find the previous and next al
    alnp = [altot(find(pp.al(1,i)<altot,1,'first')) pp.al(1,i) altot(find(pp.al(1,i)>altot,1,'last')) ; ...
            altot(find(pp.al(2,i)>altot,1,'last')) pp.al(2,i) altot(find(pp.al(2,i)<altot,1,'first'))];

    if size(alnp,2)==3 && d(i)>thn
    % get the point
    m = al2m(alnp(:)',sh);

    mv1 = linspaceN(m(:,5),m(:,6),ns);
    mv2 = linspaceN(m(:,1),m(:,2),ns);
    mv = linspaceN(m(:,3),m(:,4),ns);
    
    % get the 3D points
    W = GTPSW_TPS_kernelize(q_FP',mv1');
    XYZ1 = (PE*W');

    W = GTPSW_TPS_kernelize(q_FP',mv2');
    XYZ2 = (PE*W');
    
    W = GTPSW_TPS_kernelize(q_FP',mv');
    XYZ = (PE*W');

    % compute the 3D vectors
    n = (diff(XYZ,1,2));
    n1 = (cross(n,XYZ1(:,1:end-1)-XYZ(:,1:end-1)));
    n2 = (cross(XYZ2(:,1:end-1)-XYZ(:,1:end-1),n));
    
    ind = find(sqrt(sum(n.^2))<thn|sqrt(sum(n1.^2))<thn|sqrt(sum(n2.^2))<thn);
    n(:,ind) = [];
    n1(:,ind) = [];
    n2(:,ind) = [];
    
    % evaluate the angle
    betmp = zeros(1,size(n,2));
    for j=1:size(n,2)
        betmp(j) = angleVect(n1(:,j),n2(:,j),n(:,j));
    end

    % set the bending angle
    pp.be(i) = mean(betmp)/d(i);
    be(i) = pp.be(i);
    beWS(i) = mean(betmp);
    else
        pp.be(i) = 0;
            be(i) = pp.be(i);
    beWS(i) = 0;
    end
end

pp.be(isnan(pp.be)) = 0;
[pp.al , pp.be] = uncross(pp.al , pp.be);
