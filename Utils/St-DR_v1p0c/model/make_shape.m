% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function sh = make_shape(r)
%
%       sh = make_shape(r)
% 
% make a rectangular shape structure from the ration width/length
% Inputs:
% - the polygon p.
%
% Outputs
% - the shape structure sh, with fields:
%   - sh.c: corners of the shape
%   - sh.poly: polygon of the shape
%   - sh.bk, sh.c_x, sh.c_y: parameters of the arc length to 2D point function

l = 1/(2*(r+1));
w = r*l;

% corners of the shape
sh.c = [-w/2 , w/2 , w/2 , -w/2 ; -l/2 , -l/2 , l/2 , l/2];


sh.poly = [sh.c sh.c(:,1)];

d = sqrt(sum(diff(sh.poly , 1 , 2).^2));

sh.bk = [0 cumsum(d)];
sh.c_x = [diff(sh.poly(1,:))'./d' sh.c(1,:)'];
sh.c_y = [diff(sh.poly(2,:))'./d' sh.c(2,:)'];

sh.ppx = mkpp(sh.bk, sh.c_x);
sh.ppy = mkpp(sh.bk, sh.c_y);