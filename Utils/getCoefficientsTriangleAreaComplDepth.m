
% Check ./SDPExperiments/CalcEquiarealComplDepth.mw (eqn. 3) for the 23 monomial
% bases for which the coefficients are returned here:
function [cg] = getCoefficientsTriangleAreaComplDepth(G, pInfty_j, pInfty_q, pInfty_r)

    Delta = [pInfty_j pInfty_q pInfty_r];
    j = 1; q = 2; r = 3;
    
    cg(1) = Delta(j)^2*Delta(r)^2*G(3)/4 + Delta(q)^2*Delta(r)^2*G(6)/4 + Delta(j)^2*Delta(q)^2*G(1)/4 + Delta(j)^2*Delta(q)*Delta(r)*G(2)/4 + Delta(j)*Delta(q)^2*Delta(r)*G(4)/4 + Delta(j)*Delta(q)*Delta(r)^2*G(5)/4;
    cg(2) = G(4)/4;
    cg(3) = G(5)/4;
    cg(4) = -Delta(r)*G(5)/2 - Delta(j)*G(2)/2 - Delta(q)*G(4)/2;
    cg(5) = G(2)/4;
    cg(6) = G(6)/4;
    cg(7) = G(1)/4;
    cg(8) = G(3)/4;
    cg(9) = Delta(q)*Delta(r)*G(5)/2 + Delta(j)*Delta(q)*G(2)/2 + Delta(j)*Delta(r)*G(3) + Delta(q)^2*G(4)/4;
    cg(10) = Delta(j)*Delta(q)*G(1) + Delta(j)*Delta(r)*G(2)/2 + Delta(q)*Delta(r)*G(4)/2 + Delta(r)^2*G(5)/4;
    cg(11) = -Delta(q)*G(1)/2 - Delta(r)*G(2)/4;
    cg(12) = -Delta(j)*G(1)/2 - Delta(r)*G(4)/4;
    cg(13) = -Delta(q)*G(2)/4 - Delta(r)*G(3)/2;
    cg(14) = -Delta(j)*G(3)/2 - Delta(q)*G(5)/4;
    cg(15) = -Delta(j)*G(4)/4 - Delta(r)*G(6)/2;
    cg(16) = -Delta(j)*G(5)/4 - Delta(q)*G(6)/2;
    cg(17) = Delta(q)^2*G(1)/4 + Delta(r)^2*G(3)/4 + Delta(q)*Delta(r)*G(2)/4;
    cg(18) = Delta(j)*Delta(q)*G(4)/2 + Delta(j)*Delta(r)*G(5)/2 + Delta(q)*Delta(r)*G(6) + Delta(j)^2*G(2)/4;
    cg(19) = Delta(j)^2*G(3)/4 + Delta(q)^2*G(6)/4 + Delta(j)*Delta(q)*G(5)/4;
    cg(20) = Delta(j)^2*G(1)/4 + Delta(r)^2*G(6)/4 + Delta(j)*Delta(r)*G(4)/4;
    cg(21) = -Delta(j)*Delta(q)*Delta(r)*G(2)/2 - Delta(j)*Delta(q)^2*G(1)/2 - Delta(j)*Delta(r)^2*G(3)/2 - Delta(q)^2*Delta(r)*G(4)/4 - Delta(q)*Delta(r)^2*G(5)/4;
    cg(22) = -Delta(j)*Delta(q)*Delta(r)*G(4)/2 - Delta(j)^2*Delta(q)*G(1)/2 - Delta(j)^2*Delta(r)*G(2)/4 - Delta(j)*Delta(r)^2*G(5)/4 - Delta(q)*Delta(r)^2*G(6)/2;
    cg(23) = -Delta(j)*Delta(q)*Delta(r)*G(5)/2 - Delta(j)^2*Delta(q)*G(2)/4 - Delta(j)^2*Delta(r)*G(3)/2 - Delta(j)*Delta(q)^2*G(4)/4 - Delta(q)^2*Delta(r)*G(6)/2;

end