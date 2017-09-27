function [waveSpeed,pearson] = GetWaveSpeed(time, maxDistance, isPlot, iBeg, iEnd)
% GETWAVESPEED Linear fit to the spread of the wave of infection
% Fits from (time to maximal spread)/6 to (time to maximal spread)*3/4

if (nargin < 4 && nargin > 1)
    maxMaxDistance = max(maxDistance);
    iBeg = find( maxDistance > maxMaxDistance*1/5, 1 );
    iEnd = find( maxDistance > maxMaxDistance*4/5, 1 );
    iPlotEnd = find(maxDistance == maxMaxDistance, 1);
end

p = polyfit(time(iBeg:iEnd), maxDistance(iBeg:iEnd), 1);
waveSpeed = p(1);

pearson = corrcoef(time(iBeg:iEnd), maxDistance(iBeg:iEnd));
pearson = pearson(1,end);

if (nargin > 2 && isPlot == 1)
    hold on;
    plot(time(1:iPlotEnd), maxDistance(1:iPlotEnd), 'k-');
    plot(time(iBeg:iEnd), polyval(p, time(iBeg:iEnd)), 'r-', 'LineWidth', 2);
    hold off;
    xlabel('Time [year]');
    ylabel('Distance to disease front [km]');
    legend('Simulation','Fit','Location','Best');
    
    fprintf('R^2: %.4f\n', pearson);
end
