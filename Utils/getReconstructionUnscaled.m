function [reconstruction, scaleList] = getReconstructionUnscaled(Data, delta, dirVecs, nFiles)
    for iF = 1:nFiles
        pts = dirVecs{iF}.*delta(iF,:);
        reconstruction{iF} = pts.';
        scale{iF} = 1;
    end
    scaleList = cell2mat(scale);
end