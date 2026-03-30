function [reconstruction, scaleList] = getReconstruction(Data, delta, dirVecs, nFiles)
    for iF = 1:nFiles
        pts = dirVecs{iF}.*delta(iF,:);
        [regParams] = absor(pts, Data.Pgth(iF).P,'doScale',true,'doTrans', false); 
        reconstruction{iF} = regParams.s.*pts.';
        scale{iF} = regParams.s;
    end
    scaleList = cell2mat(scale);
end