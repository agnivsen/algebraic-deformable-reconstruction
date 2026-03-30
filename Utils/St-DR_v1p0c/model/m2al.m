% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function al = m2al(m , sh)


option= optimset('display' , 'off', 'LargeScale', 'on', 'MaxIter', 20, 'MaxFunEvals', 100000);
n_samples = 100;
al_ref = linspace(0 , 1 , n_samples);
m_ref = al2m(al_ref , sh);

for i = 1:size(m,2)
    [nth , ind] = min(sum((m_ref - repmat(m(:,i),1,n_samples)).^2));
    al(i) = lsqnonlin(@(al_i)d_sh2m(al_i , sh , m(:,i)) , al_ref(ind) , [] , [] , option);
end

return

function e = d_sh2m(al, sh , m)
e = m - al2m(al , sh);
return