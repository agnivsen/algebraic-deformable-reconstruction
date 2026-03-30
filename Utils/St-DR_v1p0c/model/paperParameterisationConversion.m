% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function pp_out = paperParameterisationConversion(pp_in,type_out,n_extra)


if ~exist('n_extra','var')
    n_extra = 0;
end

% conversion to an al-be parameterisation
switch pp_in.type
    case 'al_be'
        pp_al_be = pp_in;
    case 'valid_al_be'
        pp_al_be = valid2al(pp_in);
    case 'graphics'
        pp_al_be = graphics2al(pp_in);
end

% order the parameter
[pp_al_be,rg] = oabtr(pp_al_be);

% add extra rules if needed
if n_extra > 0
    pp_al_be = interpRules(pp_al_be,n_extra,rg);
end

% conversion to the output type
switch type_out
    case 'al_be'
        pp_out = pp_al_be;
    case 'valid_al_be'
        pp_out = al2valid(pp_al_be);
    case 'graphics'
        pp_out = al2graphics(pp_al_be,rg);
end
pp_out.r = pp_in.r;


%%
function pp_out = valid2al(pp_in)
ald = pp_in.ald;
ind = pp_in.ind;
al_one = ald(1);
ald = ald(2:end);
al = atan(ald)/pi+1/2;
al = al./sum(al);
al = [al_one al(1:end-1)];
al = mod(cumsum(al),1);
al(ind) = al;
al = reshape(al,2,length(al)/2);
pp_out.be = pp_in.be;
pp_out.al = al;
pp_out.type = 'al_be';

function pp_out = al2valid(pp_in)
pp_out.be = pp_in.be;
[al,ind] = sort(pp_in.al(:)');
ald = diff([al al(1)+1]);
ald = ald./sum(ald);
ald = tan(pi*(ald-1/2));
ald = [al(1) ald];
pp_out.ald = ald;
pp_out.ind = ind;
pp_out.type = 'valid_al_be';

%%
function pp_out = graphics2al(pp_in)

ind = pp_in.ind_rules~=0;
al = pp_in.al(ind);
be = pp_in.be(ind);
ind_rules = pp_in.ind_rules(ind);
[nth,ind] = sort(ind_rules);
n = length(al)/2;
al = reshape(al(ind),2,n);
be = mean(reshape(be(ind),2,n));
pp_out.al = al;
pp_out.be = be;
pp_out.type = 'al_be';


function pp_out = al2graphics(pp_in,rg)
n = length(pp_in.be);

% compute the number of rules per region
rgrules = zeros(max(rg),n); % a binary matrix (#region x #rules): 1 if the rule belongs to the region
rgrules(sub2ind(size(rgrules),rg,1:n)) = 1; 
nrpr = sum(rgrules,2)'; % number of rules per region

% get the end and link points:
al = pp_in.al; % arc length of the rules
be = repmat(pp_in.be,2,1); % be of the rules
rules = repmat(1:n,2,1); % rules indices
regions = repmat(rg,2,1); % region indices
sgnrules = [ones(1,n) ; -ones(1,n)]; % sign of rule points: 1 for the first point, 0 for the second one

deltaRules = 1./repmat(sum(rgrules .* repmat(nrpr',1,n),1),2,1); % distance between rules (normalisation of region length)

% order the intersection points
[al,ind] = sort(al(:)');
be = be(ind);
rules = rules(ind);
regions = regions(ind);
sgnrules = sgnrules(ind);
deltaRules = deltaRules(ind);
deltaRules2 = deltaRules;

% find the end points, the link points and the first/last points
ind_end = find(diff(rules)==0)+.5;
delta_end = deltaRules(ind_end+.5)/2;
deltaRules2(ind_end-.5) = deltaRules(ind_end-.5)/2;
sgn_end = sgnrules(ind_end+.5);

ind_rg = find(diff(regions)~=0)+.5;
delta_rg = deltaRules(ind_rg+.5)/2;
deltaRules2(ind_rg-.5) = deltaRules(ind_rg-.5)/2;
sgn_rg = sgnrules(ind_rg+.5);

ind_fl = [.5 2*n+.5];
delta_fl = deltaRules([1 end])/2;
deltaRules2(end) = deltaRules(end)/2;
sgn_fl = [1 -1];

deltaRules = deltaRules2;

% identification of the points
% 0: rule -- 1: first last (fl) -- 2: end -- 3: link
id_rules = zeros(size(al));
id_fl = [1 1];
id_end = 2*ones(size(ind_end));
id_rg = 3*ones(size(ind_rg));

% compute the al and the be for the limit points
alpchip = pchip(0:2*n+1,[al(end)-1 al al(1)+1]);
bespline = spline(0:2*n+1,[be(end) be be(1)]);

al_rg = ppval(alpchip,ind_rg);
be_rg = ppval(bespline,ind_rg);

al_end = ppval(alpchip,ind_end);
be_end = ppval(bespline,ind_end);

al_fl = ppval(alpchip,ind_fl);
be_fl = ppval(bespline,ind_fl);

% merge all data (rule, end, link, first/last) in vectors
al_tot = [al al_rg al_end al_fl];
be_tot = [be be_rg be_end be_fl];
delta_tot = [deltaRules delta_rg delta_end delta_fl];
id_tot = [id_rules id_rg id_end id_fl];
sgn_tot = [sgnrules sgn_rg sgn_end sgn_fl];
rules_tot = [rules zeros(size(id_rg)) zeros(size(id_end)) zeros(size(id_fl))];

% order the vectors against points arc lengthes
[al_tot,ind] = sort(al_tot);
be_tot = be_tot(ind);
delta_tot = delta_tot(ind);
id_tot = id_tot(ind);
sgn_tot = sgn_tot(ind);
rules_tot = rules_tot(ind);

pos = [0 cumsum(delta_tot(1:end-1).*sgn_tot(1:end-1))]; % point positions


pp_out.al = al_tot;
pp_out.be = be_tot;
pp_out.pos = pos;
pp_out.ind_rules = rules_tot;
pp_out.id = id_tot;
pp_out.type = 'graphics';


%%
function pp_out = interpRules(pp_in,n_extra,rg_in)

fr_al_in = pp_in.al(:,1);

n = length(pp_in.be);
% al interpolation
[altmp,ind] = sort(pp_in.al(:)');
indr = 1:2*n;
indr(ind) = indr;
indr = reshape(indr,2,n);

d = 0:1/(2*n_extra+1):.5;
d = [-d(end:-1:2) d];
d = [d;-d];
indrtot = repmat(indr,[1 1 2*n_extra+1]) + repmat(reshape(d,[2 1 2*n_extra+1]),1,n);
al = pchip(-1:2*n+2,[altmp([end-1 end])-1 altmp altmp([1 2])+1],indrtot);
pp_out.al = reshape(al,2,n*(2*n_extra+1));

% if the first region has changed due to the new rules, the following code changes the sign of bending angles
pp_out.be = zeros(1,size(pp_out.al,2));
[pp_out_tmp,rg] = oabtr(pp_out);


%%%%%%%%%%
ind1 = find(mod(fr_al_in(1),1)==mod(pp_out_tmp.al(1,:),1) & mod(fr_al_in(2),1)==mod(pp_out_tmp.al(2,:),1) |...
            mod(fr_al_in(2),1)==mod(pp_out_tmp.al(1,:),1) & mod(fr_al_in(1),1)==mod(pp_out_tmp.al(2,:),1));

if rg(ind1)~=1
    pp_in.be(rg_in==1) = -pp_in.be(rg_in==1);
    ind2 = find(mod(pp_out_tmp.al(1,n_extra+1),1)==mod(pp_in.al(1,:),1) & mod(pp_out_tmp.al(2,n_extra+1),1)==mod(pp_in.al(2,:),1) |...
                mod(pp_out_tmp.al(2,n_extra+1),1)==mod(pp_in.al(1,:),1) & mod(pp_out_tmp.al(1,n_extra+1),1)==mod(pp_in.al(2,:),1));
    pp_in.be(rg_in==rg_in(ind2)) = -pp_in.be(rg_in==rg_in(ind2));
end

% be interpolation
betmp = [pp_in.be ; pp_in.be];
betmp = vect(betmp(ind))';
be = pchip(0:2*n+1,[betmp(end) betmp betmp(1)],indrtot);

pp_out.be = mean(reshape(be,2,n*(2*n_extra+1)));
pp_out = oabtr(pp_out);


