% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function teta_OK = find_angles(R)

if norm(R'*R)-1 > 1e-8
    error('R is not a rotation matrix');
else
    teta(2) = -asin(R(3,1));

    teta(3) = atan(R(3,2) / R(3,3));
    teta(1) = atan(R(2,1)/R(1,1));

    teta = [teta(:)' ; [pi + teta(1),  pi - teta(2),  pi + teta(3)]];
    
    for i = 1:2
        teta_p(1) =teta(i,1);
        for j = 1:2
            teta_p(2) =teta(j,2);
            for k = 1:2
                teta_p(3) =teta(k,3);
                Rp = rotation(teta_p);
                if norm(Rp-R) < 1e-8
                    teta_OK = teta_p;
                end
            end
        end
    end
end
