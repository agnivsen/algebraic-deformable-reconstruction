% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function paper = randPaper(r)

% it could be better to draw all the random values at the same time

n_guidingMax = 5; % the maximal number of guiding rulings per regions
n_interpolation = 4; % the interpolation ruling number
if ~exist('r','var')
    pp.r = randInt(.5,1.5);
else
    pp.r = r;
end
alpha_min = pi/(12*n_interpolation);
alpha_max = pi/(12*n_interpolation);

% choose the number of regions on the paper
n_regions = [1 3 4 5];
n_regions = n_regions(randInteger(1,4));

% choose the number of guiding rulings per regions
n_guiding = randInteger(1,n_guidingMax,[1,n_regions]);
n_rulings = sum(n_guiding);

% draw the ruling arc lengths -- not the good method, first draw the limit with constraints (such that
% the regions are not only along the boundary) and then draw the arc lengths

% draw the hot points and region arc lenghts
sh = make_shape(pp.r);
switch n_regions
    case 1
        hot_points = sort(rand(1,2));
        
        pp.al = [sort(randInt(hot_points(1),hot_points(2),1,n_rulings)) ; ...
                 sort(randInt(hot_points(2),hot_points(1)+1,1,n_rulings),'descend')];
    case 3
        ind_border = randperm(4);
        hot_points = zeros(1,3);
        for i=1:3
            hot_points(i) = randInt(sh.bk(ind_border(i)),sh.bk(ind_border(i)+1));
        end
        hot_points = sort(hot_points);
        hot_points = [hot_points hot_points(1)+1];
        pp.al = [];
        for i=1:3
            altmp = sort(randInt(hot_points(i),hot_points(i+1),1,2*n_guiding(i)));
            pp.al = [pp.al [altmp(1:n_guiding(i)) ; altmp(end:-1:n_guiding(i)+1)]];
        end
    case 4
        hot_points = zeros(1,4);
        for i=1:4
            hot_points(i) = randInt(sh.bk(i),sh.bk(i+1));
        end
        hot_points = sort(hot_points);
        hot_points = [hot_points hot_points(1)+1];
        pp.al = [];
        for i=1:4
            altmp = sort(randInt(hot_points(i),hot_points(i+1),1,2*n_guiding(i)));
            pp.al = [pp.al [altmp(1:n_guiding(i)) ; altmp(end:-1:n_guiding(i)+1)]];
        end
    case 5
        dbl = randInteger(1,2);
        hot_points = zeros(1,6);
        % for dbl1
        hot_points(1:2) = randInt(sh.bk(dbl),sh.bk(dbl+1),1,2);
        hot_points(3:4) = randInt(sh.bk(dbl+2),sh.bk(dbl+3),1,2);
        hot_points(5) = randInt(sh.bk(dbl+1),sh.bk(dbl+2));
        if dbl==1
            hot_points(6) = randInt(sh.bk(4),sh.bk(5));
        else
            hot_points(6) = randInt(sh.bk(1),sh.bk(2));
        end
        hot_points = sort(hot_points);
        
        if dbl==1
            hot_points = [hot_points hot_points(1)+1];
        else
            hot_points = [hot_points(2:end) hot_points(1:2)+1];
        end
        pp.al = [sort(randInt(hot_points(1),hot_points(2),1,n_guiding(1))) ; sort(randInt(hot_points(4),hot_points(5),1,n_guiding(1)),'descend')];
        ind = [2 3 5 6 ; 3 4 6 7];
        for i=1:4
            altmp = sort(randInt(hot_points(ind(1,i)),hot_points(ind(2,i)),1,2*n_guiding(i+1)));
            pp.al = [pp.al [altmp(1:n_guiding(i+1)) ; altmp(end:-1:n_guiding(i+1)+1)]];
        end
end

% draw the bending angles
pp.be = randInt(alpha_min, alpha_max,1,n_rulings);
pp.type = 'al_be';
pp = paperParameterisationConversion(pp,'al_be',n_interpolation);

D = newPaperComputeD(pp,sh);
pp.be(D~=0) = pp.be(D~=0)./D(D~=0);
pp.be(D==0) = 0;
paper = newPaperFast(pp);
