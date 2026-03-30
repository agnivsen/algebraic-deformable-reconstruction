% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function p = plotMesh(meshIn , varargin)
%
% plotMesh(meshIn , varargin)
%
% plot the mesh meshIn with 'varargin' plot options
if ishold
    c=1;
else
    c=0;
end

if size(meshIn , 3) == 3
    p1 = plot3(meshIn(:,:,1) , meshIn(:,:,2) , meshIn(:,:,3) , varargin{:});
    hold on
    p2 = plot3(meshIn(:,:,1)' , meshIn(:,:,2)' , meshIn(:,:,3)' , varargin{:});
elseif size(meshIn , 3) == 2    
    p1 = plot(meshIn(:,:,1),meshIn(:,:,2) , varargin{:});
    hold on
    p2 = plot(meshIn(:,:,1)',meshIn(:,:,2)', varargin{:});   
else
    error('mesh non conforme')
end

if c==0
    hold off
end

if nargout>0
p = [p1 p2];
end