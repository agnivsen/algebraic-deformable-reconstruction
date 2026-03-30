function [Data, nng, template, K] =  loadHulkData(nK, m)
    data = load('./Data/data_hulk.mat');

    flatIndex = 8; % the 8th frame of this dataset is actually flat (serves as template)

    if ~exist('m', 'var')
        m = 1000;
    end
    
    mActual = size(data.Pgth(1).P,2);

    if(mActual < m)
        indices = 1:mActual;
    else 
        indices = randperm(mActual);
        indices = sort(indices(1:m));
    end

    K = load('./Data/intrinsics_hulk.txt');
    
    v = [];
    for iN = 1:size(data.Pgth,2)
        Data.p(iN).p = data.p(iN).p(:,indices);
        Data.Pgth(iN).P = data.Pgth(iN).P(:,indices);
        v = [v; ones(1,numel(indices))];
    end
    Data.v = v;
    
    template = data.Pgth(flatIndex).P(:,indices).';

    template = template./1000; % doing stuff in metres

    nFiles = size(Data(1).p,2);

    for ii = 1:nFiles
        Data.Pgth(ii).P = Data.Pgth(ii).P./1000;
    end
    
    [nng] = getNeighborhoodDuplicated(Data,nK);

    fprintf('Running <strong>HULK</strong> data with: nPts = %d\n\n', numel(indices));
end