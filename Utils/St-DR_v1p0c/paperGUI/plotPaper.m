% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function plotPaper



global rulings handlesG pp figure_3D paper randomPoints_2D randomPoints_3D mesh3D mesh2D T
fbe = 20; % coefficient to apply on the angle to get a good deformation range

if size(rulings,2)>1 % it is a bug in the model: it does not work whith only one ruling
    pp.al = rulings(1:2,:);
    pp.be = rulings(3,:)*fbe;
    pp.type = 'al_be';
    n_extra = str2double(get(handlesG.model_editExtra,'string'));
    pp = paperParameterisationConversion(pp,'al_be', n_extra);

    paper = newPaperFast(pp);

    nv = str2double(get(handlesG.display_editMesh,'string'));
    [M,m] = paperMesh(paper,nv);
    %[M,T] = meshRigidTransformation(M,cat(3,m,zeros(size(m,1),size(m,2))));
    %This was the Perriollat's original code.
    %The new one (no autoscaling) is replaced below  
    [M,T] = meshRigidTransformation(M,cat(3,m,zeros(size(m,1),size(m,2))),0);
    % This does not scale the data automatically.
    % Modified by Jae-Hak Kim, 2011-06-09 by the help from M. Perriollat.
    
    mesh3D = M;
    mesh2D = m;
    
    figure(figure_3D)
    v = view;
    clf
    hold on
    % display mesh
    if get(handlesG.display_radiobuttonMesh,'value')
        plotMesh(M,'k')        
        
        hold on
%         load paperRef_mesh
%         plotMesh(mesh3D,'r')
    end
    
    % display rulings
    if get(handlesG.display_radiobuttonRulings,'value')
        n = length(pp.be);        
        R = T*[uv2M([paper.uv_rules(1:2,:) paper.uv_rules(3:4,:)],paper);ones(1,2*n)];
        plot3([R(1,1:n);R(1,n+1:end)],[R(2,1:n);R(2,n+1:end)],[R(3,1:n);R(3,n+1:end)],'m','linewidth',3)
        
    end
    
    if get(handlesG.display_radiobuttonPoints,'value') && ~isempty(randomPoints_2D)
        rp(1,:) = (paper.sh.c(1,2)-paper.sh.c(1,1))*randomPoints_2D(1,:) + paper.sh.c(1,1);
        rp(2,:) = (paper.sh.c(2,3)-paper.sh.c(2,2))*randomPoints_2D(2,:) + paper.sh.c(2,2);
        randomPoints_3D = T*[uv2M(rp,paper);ones(1,size(rp,2))];
        plot3(randomPoints_3D(1,:),randomPoints_3D(2,:),randomPoints_3D(3,:),'b.','markersize',10)
    end
    
    axis equal
    view(v);
end