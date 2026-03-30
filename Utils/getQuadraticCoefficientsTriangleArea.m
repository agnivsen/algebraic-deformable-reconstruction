function [G] = getQuadraticCoefficientsTriangleArea(Pj, Pq, Pr)

    if( ( ~isa(Pj,'sdpvar') ) && ( ~isa(Pq,'sdpvar') ) && ( ~isa(Pr,'sdpvar') ) )
        assert(abs(norm(Pj) - 1) < 10^(-10),'getQuadraticCoefficientsTriangleArea: inputs (first argument) must be an unit normal' );
        assert(abs(norm(Pq) - 1) < 10^(-10),'getQuadraticCoefficientsTriangleArea: inputs (second argument) must be an unit normal' );
        assert(abs(norm(Pr) - 1) < 10^(-10),'getQuadraticCoefficientsTriangleArea: inputs (third argument) must be an unit normal' );

        G = zeros(1,6);
    end

    xj = Pj(1);   yj = Pj(2);   zj = Pj(3);
    xq = Pq(1); yq = Pq(2); zq = Pq(3);
    xr = Pr(1);  yr = Pr(2);   zr = Pr(3);

% % %     G(1) = ((yq^2 + zq^2)*xj^2) - (2*( (yq * xq * yj) + (zq * xq * zj) ) * xj) + ((zq^2 + xq^2) * yj^2) - ... 
% % %                                                                                                             (2 * zq * yq * zj * yj) + ((yq^2 + xq^2) * zj^2);
% % % 
% % %     G(2) = (2 * ( ( ((xr  * yq) + (xq * yr)) * yj) + (( (xr * zq) + (zr * xq)) * zj) ) * xj) -... 
% % %                             (2 * ((zr * zq) + (xq * xr)) * yj^2) + (2 * ( (zr * yq) + (zq * yr) ) * zj * yj) - (2 * ((yq * yr) + (xq * xr)) * zj^2) - (2 * ( (yq * yr) + (zr * zq) ) * xj^2);
% % % 
% % %     G(3) = ((yr^2 + zr^2) * xj^2) - (2 * ((yr * xr * yj) + (zr * xr * zj)) * xj) + ((xr^2 + zr^2) * yj^2) - (2 * yj * zj * yr * zr) ...
% % %                                                                             + (zj^2 * (xr^2 + yr^2));
% % % 
% % %     G(4) = (2 * ( ( (yq * yr) + (zr * zq) ) * xq) - (xr * (yq^2 + (zq)^2)) * xj) + (2 * ((zr * zq * yq) + (xr * yq * xq) - (zq^2 * yr) - (xq^2 * yr) ) * yj) ...
% % %                                                     - (2 * ( (yq^2 * zr) - (zq * yr * yq) - (xr * zq * xq) + (xq^2 * zr) ) * zj);
% % % 
% % %     G(5) = (2 * (  (xr * ( (yq * yr) + (zr * zq))) - ((yr^2 + zr^2) * xq) ) * xj) + (2 * ((xr * yr * xq) - ((xr^2 + zr^2) * yq) + (zq * yr * zr)) * yj) ...
% % %                                                                     + (2 * ((xr * zr * xq) + (yr * zr * yq) - (zq * (xr^2 + yr^2))) * zj);
% % % 
% % %     G(6) = ((yr^2 + zr^2) * xq^2) - (2 * xr * ( (yq * yr) + (zr * zq) ) * xq) + ((xr^2 + zr^2) * yq^2) - (2 * zq * yr * zr * yq) + (zq^2 * (xr^2 + yr^2));
    
    
    G(1) = ((yq) ^ 2 + (zq) ^ 2) * (xj) ^ 2 - 2 * (yq * xq * yj + zq * xq * zj) * xj +... 
                                ((zq) ^ 2 + (xq) ^ 2) * (yj) ^ 2 - 2 * zq * yq * zj * yj + ((yq) ^ 2 + (xq) ^ 2) * (zj) ^ 2;
                            
    G(2) = (-2 * yq * yr - 2 * zr * zq) * (xj) ^ 2 - 2 * ((-yq * xr - xq * yr) * yj + (-zq * xr - xq * zr) * zj) * xj ...
                                + (-2 * zr * zq - 2 * xr * xq) * (yj) ^ 2 - 2 * (-yq * zr - zq * yr) * zj * yj + (-2 * yq * yr - 2 * xr * xq) * (zj) ^ 2;
                            
    G(3) = ((yr) ^ 2 + (zr) ^ 2) * (xj) ^ 2 - 2 * (yr * xr * yj + zr * xr * zj) * xj + ((xr) ^ 2 + (zr) ^ 2) * (yj) ^ 2 - 2 * yj * zj * yr * zr ...
                                                                                                        + (zj) ^ 2 * ((xr) ^ 2 + (yr) ^ 2);
                                                                                                    
    G(4) = 2 * ((yq * yr + zr * zq) * xq - xr * ((yq) ^ 2 + (zq) ^ 2)) * xj - 2 * (-zr * zq * yq - xr * yq * xq + (zq) ^ 2 * yr +... 
                                                                        (xq) ^ 2 * yr) * yj - 2 * ((yq) ^ 2 * zr - zq * yr * yq - xr * zq * xq + (xq) ^ 2 * zr) * zj;
                                                                    
    G(5) = 2 * (-((yr) ^ 2 + (zr) ^ 2) * xq - xr * (-yq * yr - zr * zq)) * xj - 2 * (-xr * yr * xq + ((xr) ^ 2 + (zr) ^ 2) * yq - zq * yr * zr) * yj ...
                                                                                                                    - 2 * (-xr * zr * xq - yr * zr * yq + zq * ((xr) ^ 2 + (yr) ^ 2)) * zj;
                                                                                                                
    G(6) = ((yr) ^ 2 + (zr) ^ 2) * (xq) ^ 2 - 2 * xr * (yq * yr + zr * zq) * xq + ((xr) ^ 2 + (zr) ^ 2) * (yq) ^ 2 - 2 * zq * yr * zr * yq ...
                                                                                                                                                                    + (zq) ^ 2 * ((xr) ^ 2 + (yr) ^ 2);
end