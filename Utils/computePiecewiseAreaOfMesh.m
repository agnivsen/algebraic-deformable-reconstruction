function [area] = computePiecewiseAreaOfMesh(vertices, nng)
    nPoints = size(nng,1);
    nNeighbors = size(nng,2);

    assert(nNeighbors > 1,'computePiecewiseAreaOfMesh: area cannot be computed for less than 2 neighbors');

    assert(((size(vertices,2) == 2) || (size(vertices,2) == 3)) ,...
                        'computePiecewiseAreaOfMesh: area can be computed only for 2 or 3 dimensional vertices');

    for iPoint = 1:nPoints
        for iNeighbor = 1:(nNeighbors-1)
            P1 = vertices(iPoint,:);
            P2 = vertices(nng(iPoint, iNeighbor),:);
            P3 = vertices(nng(iPoint, (iNeighbor + 1),:));

            triangleIndex = iNeighbor;

            area{iPoint, triangleIndex} = areaOfTriangle(P1, P2, P3);
        end
    end
end