function [rmseList] = plotResultsNRSfM(reconstruction, data)

    M = size(data.Pgth,2);

    rmseList = []; increment = 150;

    xPos = 0;

    for ii = 1:M
        f = figure(ii);
        f.Position = [xPos 0 300 300];
        xPos = xPos + increment;

        [rms] = plotResults(reconstruction{ii}, data.Pgth(ii).P.');
        rmseList = [rmseList rms];
        
        pause(0.01);
    end


end