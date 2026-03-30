% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function [x_new,y_new] = enforceDecreasing(x,y)

% keep only the points that satisfy the decreasing condition
xy = sortrows([x(:) y(:)],[1 -2]);

[nth,ind_x] = sort(xy(:,1));
[nth,ind_y] = sort(xy(:,2),'descend');
ind = ind_x(ind_x==ind_y);
x_filtered = xy(ind,1);
y_filtered = xy(ind,2);

x_new = [];
y_new = [];

if length(x_filtered)>1
    % interpolate to fill the holes and filter to avoid multiple neighboor (8-connexity)
    xy_diff = diff([x_filtered y_filtered]);
    [i,j] = find(~(xy_diff==0 | abs(xy_diff)==1));

    for k=1:length(i)
        if j(k)==1
            x_new_tmp = x_filtered(i(k)):x_filtered(i(k)+1);
            y_new_tmp = linspace(y_filtered(i(k)),y_filtered(i(k)+1),length(x_new_tmp));
        else
            y_new_tmp = y_filtered(i(k)):-1:y_filtered(i(k)+1);
            x_new_tmp = linspace(x_filtered(i(k)),x_filtered(i(k)+1),length(y_new_tmp));
        end
        x_new = [x_new x_new_tmp(2:end-1)];
        y_new = [y_new y_new_tmp(2:end-1)];
    end

    xy_new = sortrows([x_new' y_new' ; x_filtered y_filtered],[1 -2]);
    xy_new = xy_new - repmat(mean(xy_new),length(xy_new),1) + repmat(mean(xy),length(xy_new),1);
    x_new = xy_new(:,1);
    y_new = xy_new(:,2);
end



