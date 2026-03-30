function [area] = areaOfTriangle(P1, P2, P3)
    d1 = EuclideanDistance(P1, P2);
    d2 = EuclideanDistance(P2, P3);
    d3 = EuclideanDistance(P3, P1);

    s = (d1 + d2 + d3)./2;

    area = sqrt(s .* (s - d1) .* (s - d2) .* (s - d3));
end