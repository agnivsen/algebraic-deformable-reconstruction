% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function [M,m,ind] = paperMesh(paper,varargin)

if isempty(varargin)
    nv = 10;
elseif length(varargin)==1
    nv = varargin{1};
else
    a = varargin{1};
    b = varargin{2};
end

if ~exist('a','var')
    a = round(nv*sqrt(paper.ppnew.r));
    b = round(nv/sqrt(paper.ppnew.r));
end

x = linspace(paper.sh.c(1,1),paper.sh.c(1,2),a);
y = linspace(paper.sh.c(2,1),paper.sh.c(2,3),b);

[X,Y] = meshgrid(x,y);

m = cat(3,X,Y);


[M,ind] = uv2M([X(:) Y(:)]',paper);
M = reshape(M',[b,a,3]);
ind = reshape(ind,[b,a]);
