% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function imLabel3 = cutAndGuess(imLabel)
% working implementation, but not efficient at all

imLabel2 = cutComponents(imLabel);
% cut
imLabel2 = bwlabel(imLabel2~=0);
n_curves = max(imLabel2(:));
n_samples = size(imLabel,1);
imTest = zeros(n_samples);
for k=1:n_curves
    [i,j] = find(imLabel2==k);
    p = [min(j) 1 1 min(j) max(j)+1 max(j)+1 min(i)-1 max(i) max(i) min(j) ;...
         max(i) max(i) min(j) min(j) max(j)+1 min(i)-1 min(i)-1 max(i) n_samples n_samples];
    imTest = imTest + roipoly(imTest,p(1,:)-1,p(2,:));
end

% get the block that can contain the new curves
imTest = tril(imTest==0,-1);

% replace each block by its diagonal
imTest2 = zeros(n_samples);
imTest = bwlabel(imTest);
n_blocks = max(imTest(:));
for k=1:n_blocks
    [i,j] = find(imTest==k);
    n = max(max(j)-min(j),max(i)-min(i))+1;
    is = round(linspace(min(i),max(i),n));
    js = round(linspace(max(j),min(j),n));
    imTest2(sub2ind([n_samples n_samples],is,js)) = 1;
end

imTest2 = bwlabel(tril((imTest2+imLabel2)~=0,-1));
imLabel3 = zeros(n_samples);
for k=1:n_curves
    [i,j] = find(imLabel2==k,1);
    [i,j] = find(imTest2==imTest2(i,j));
    imLabel3(sub2ind([n_samples n_samples],i,j)) = 1;
end
imLabel3 = bwlabel(imLabel3==1);
imLabel3 = cutComponents(imLabel3);
