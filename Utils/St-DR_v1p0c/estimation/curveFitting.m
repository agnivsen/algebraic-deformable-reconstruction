% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function [CP,alnew,benew] = curveFitting(X,Y,Z,th)

flag_figure=1; % flag to display figures
kz = .5; % weight on Z

Xold = X;
Yold = Y;
Zold = Z;

% filtering the data
n = 5;
X = medfilt1(X,n);
Y = medfilt1(Y,n);
Z = medfilt1(Z,n);

X0 = X(1);
Y0 = Y(1);
Z0 = Z(1);


[b,a] = butter(3,.1,'low');
X = filter(b,a,X-X0) + X0;
Y = filter(b,a,Y-Y0) + Y0;
Z = filter(b,a,Z-Z0) + Z0;

% normalisation of the data -- X and Y represent the same thing
X2 = X/std([X Y]);
Y2 = Y/std([X Y]);
Z2 = kz*Z/std(Z);

if ~exist('th','var')
    th = .01;
end

l = length(X);
lf = 10*l; % the sampling of the spline
f = 10;

indCP = [1 l];


CPb = [X2([1 end]) ; Y2([1 end]) ; Z2([1 end])];
CPin = [];
c=1;
if ~flag_figure
figure
plot3(X2,Y2,Z2,'b','linewidth',3)
hold on
text(-4.9,-7,17.8,'arc lengths','fontsize',30)
zlabel('bending angle','fontsize',30)
set(gca,'xticklabel','','yticklabel','','zticklabel','')
grid on
end

list = [2,5,12];
while c
    % defintion of the spline
    CP = [X2(indCP) ; Y2(indCP) ; Z2(indCP)];
    t = linspace(0,1,size(CP,2)); % the knot vector
    fx = pchip(t,CP(1,:));
    fy = pchip(t,CP(2,:));
    fz = pchip(t,CP(3,:));

    ts = linspace(0,1,lf);
    xs = ppval(fx,ts);
    ys = ppval(fy,ts);
    zs = ppval(fz,ts);

    ex = (repmat(X2,lf,1) - repmat(xs',1,l)).^2;
    ey = (repmat(Y2,lf,1) - repmat(ys',1,l)).^2;
    ez = (repmat(Z2,lf,1) - repmat(zs',1,l)).^2;
    dst = ex+ey+ez;

    for i=1:length(indCP)-1
        ind = ts<=t(i) | ts>=t(i+1);
        dst(ind,indCP(i):indCP(i+1)) = nan;
    end
    [d,ind] = min(dst);
    if ~flag_figure
    p(1) = plot3(xs,ys,zs,'r','linewidth',5);
    p(2) = plot3(CP(1,:),CP(2,:),CP(3,:),'r.','markersize',50);
    axis tight
    delete(p)
    end
    % end or control point insertion where the error is maximum
    [sc,indSC] = max(d);
    if sc<th^2
        c=0; % end of the optimisation process
    else % add a new control point        
        indCP = sort([indCP indSC]);
    end
end

fprintf('number of control points=%d\n',size(CP,2))

CP(1:2,:) = CP(1:2,:)*std([X Y]);
CP(3,:) = CP(3,:)*std(Z);

t = linspace(0,1,size(CP,2));
fx = pchip(t,CP(1,:));
fy = pchip(t,CP(2,:));
fz = pchip(t,CP(3,:));

t = linspace(0,1,lf);
if ~flag_figure
    figure
    plot3(X,Y,Z)
    hold on
    plot3(X,Y,Z,'b.','markersize',10)
    plot3(Xold,Yold,Zold,'k')
    plot3(ppval(fx,t),ppval(fy,t),ppval(fz,t),'r')
    plot3(CP(1,:),CP(2,:),CP(3,:),'r.','markersize',10)
end
t = linspace(0,1,size(X,2));
alnew = [ppval(fx,t);ppval(fy,t)];
benew = ppval(fz,t);
