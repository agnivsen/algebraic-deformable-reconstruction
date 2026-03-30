% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function paper = newPaperFast(pp,st)

% pp.r: ratio lenght/height
% pp.al: rule arc lengths
% pp.be: bending
% st: show time flag
% paper.test: giving uv, return the region
% paper.RT: giving the region and uv, compute the 3D points

threshold; % define thresholds for geometric properties

if ~exist('st','var')
    st = 0;
elseif st == 1
    fprintf('\n#########\n');tic
end

% 1. check non-crossing
if ~isempty(get_intersections(pp.al))
    error('some rules intersect each other')
end

n = length(pp.be);
sh = make_shape(pp.r);
% al_corners = [sh.bk sh.bk+1];

if st
    fprintf('parameter checking time = %f\n',toc);tic;
end
% 2. order and bring together regions
[ppnew,rg] = oabtr(pp);

uv_rules = reshape(al2m(ppnew.al(:)',sh),4,n);
test = eye(n);
RT = repmat([eye(3) zeros(3,1)],n+1,1);
rulesDependencies = zeros(n);

% compute the rule points position
[altmp,ind] = sort(ppnew.al(:)');
indr = 1:2*n;
indr(ind) = indr;
indr = reshape(indr,2,n);

if st
    fprintf('organising parameter time = %f\n',toc);tic;
end

d = newPaperComputeD(pp,sh);

for i=1:n
    % 3. compute test matrix
    ind1 = i+find(ppnew.al(1,i+1:end)>=ppnew.al(1,i),1,'first');
    ind2 = i+find(ppnew.al(2,i+1:end)<=ppnew.al(2,i));
    [nth,ind] = max(ppnew.al(2,ind2));
    next = [ind1 ind2(ind)];

    if length(next)>1
        if diff(next) == 0
            next = next(1);
        else
            ind = find(ppnew.al(1,:)>=ppnew.al(2,next(1)) & ppnew.al(2,:)<=ppnew.al(1,next(2)),1,'first');
            while ~isempty(ind)
                next = [next ind];
                ind = find(ppnew.al(1,:)>=ppnew.al(2,next(end)) & ppnew.al(2,:)<=ppnew.al(1,next(2)),1,'first');
            end
        end
    end
    
    
    test(i,next) = -1;
    
    % 4. form transformations
    agl = d(i)*ppnew.be(i);
    m = al2m(ppnew.al(:,i)',sh);
    
    indRT = find(ppnew.al(1,:)>=ppnew.al(1,i) & ppnew.al(2,:)<=ppnew.al(2,i));
    
    rulesDependencies(i,indRT) = 1;
    RT_i = compute_RT2(m,agl);
    ind = linspaceN(3*(indRT)+1,3*(indRT+1),3);
    RT(ind(:),:) = RT(ind(:),:) * [RT_i;zeros(1,3) 1];
end

if st
    fprintf('model computation time = %f\n',toc);tic;
end


test = [-1 zeros(1,n-1) ; test];


RT = reshape(permute(reshape(RT,[3,n+1,4]),[1 3 2]),[3,4*(n+1)]);
RTmp(diagPattern(n+1,3,4)==1) = RT;
RT = reshape(RTmp,3*(n+1),4*(n+1));

% compute matrices for the region tests
v = uv_rules(3:4,:) - uv_rules(1:2,:);
A = [v(2,:) ; -v(1,:)];
b = -sum(uv_rules(1:2,:).*A);
% % 5. set up rigid transformation
% T = 
% RotationUsingQuat

% save paper
paper.ppnew = ppnew;
paper.sh = sh;
paper.test.mask = test./repmat(sum(abs(test),2),1,n);
paper.test.A = A';
paper.test.b = b';
paper.RT = sparse(RT);
paper.uv_rules = uv_rules;
paper.rulesDependencies = rulesDependencies;
paper.rg = rg;

if st
    fprintf('model saving time = %f\n\n',toc);tic;
end

function RT = compute_RT(m,be)
% distance computation
Ar = polyarea(m(1,:) , m(2,:));
if size(m,2) == 3
    l = sqrt(sum((m(:,1)-m(:,3)).^2))/2;
else
    l = (sqrt(sum((m(:,1)-m(:,4)).^2)) + sqrt(sum((m(:,2)-m(:,3)).^2)))/2;
end
d = Ar/l;

R = rot_vect([(m(:,end) - m(:,1));0],be*d);
T = -R*[m(:,1);0] + [m(:,1);0];
RT = [R T];

function RT = compute_RT2(m,agl)
% distance computation

R = rot_vect([(m(:,2) - m(:,1));0],agl);
T = -R*[m(:,1);0] + [m(:,1);0];
RT = [R T];

