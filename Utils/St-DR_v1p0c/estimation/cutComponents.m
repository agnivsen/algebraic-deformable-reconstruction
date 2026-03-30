% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function [imLabel2,rComponents,cComponents] = cutComponents(imLabel)

nComponents = max(imLabel(:));

% get the components
rComponents = cell(1,nComponents);
cComponents = cell(1,nComponents);

for i=1:nComponents
    [r,c] = find(imLabel == i);
    [r,ind] = sort(r,'descend');
    c = c(ind);
    for j=unique(r)'
        ind = find(r==j);
        c(ind) = sort(c(ind),'ascend');
    end
    rComponents{i} = r;
    cComponents{i} = c;
end

for i=1:nComponents-1
    for j=i+1:nComponents
        go = 1;
        while go
            % get the bounds
            if ~isempty(rComponents{i})&&~isempty(rComponents{j})
                minRi = min(rComponents{i});
                maxRi = max(rComponents{i});
                minCi = min(cComponents{i});
                maxCi = max(cComponents{i});

                minRj = min(rComponents{j});
                maxRj = max(rComponents{j});
                minCj = min(cComponents{j});
                maxCj = max(cComponents{j});

                % compare
                poly_i = [maxRi size(imLabel,1) size(imLabel,1) minRi minRi minCi minCi maxRi maxRi ; minCi minCi maxRi maxRi maxCi maxCi 0 0 minCi];
                poly_j = [maxRj size(imLabel,1) size(imLabel,1) minRj minRj minCj minCj maxRj maxRj ; minCj minCj maxRj maxRj maxCj maxCj 0 0 minCj];

                ind_j = find(inpolygon(rComponents{j},cComponents{j},poly_i(1,:),poly_i(2,:)));
                ind_i = find(inpolygon(rComponents{i},cComponents{i},poly_j(1,:),poly_j(2,:)));

                % update the components
                if isempty(ind_j)&&isempty(ind_i)
                    go = 0;
                else
                    if isempty(intersect(ind_i,1)) % remove at the end of the component
                        rComponents{i}(end) = [];
                        cComponents{i}(end) = [];
                    else % at the begining
                        rComponents{i}(1) = [];
                        cComponents{i}(1) = [];
                    end
                    if isempty(intersect(ind_j,1)) % remove at the end of the component
                        rComponents{j}(end) = [];
                        cComponents{j}(end) = [];
                    else % at the begining
                        rComponents{j}(1) = [];
                        cComponents{j}(1) = [];
                    end
                end

                imLabel2 = zeros(size(imLabel));
                for k=1:nComponents
                    imLabel2(sub2ind(size(imLabel),rComponents{k},cComponents{k})) = k;
                end
                
            else
                go = 0;
            end
        end
    end
end


imLabel2 = zeros(size(imLabel));
for i=1:nComponents
    imLabel2(sub2ind(size(imLabel),rComponents{i},cComponents{i})) = i;
end