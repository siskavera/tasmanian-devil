function [r0, pearsonR0, r023, pearsonR023, eqPrev, eqPrevStd, ...
        r0Adap, pearsonR0Adap, r0Adap23, pearsonR0Adap23, eqPrevAdap, eqPrevStdAdap, ...
        r0Sigmoid, r0Sigmoid23, eqPrevSigmoid, waveSpeed, pearsonWaveSpeed, proportionReached] = ...
    Tasmanian(infectionRate, latentPeriod, contactRate, diagnoseProp, ...
    diagnosePeriod, migrationRate, xInit, yInit, K)

% Set fixed parameters
minPrev = 0.2;
minDistance = 20;  % In order to throw away runs that did not spread
minProportionReached = 0.3; % In order to throw away runs that did not spread enough

% Set return variables in case of premature extinction of disease (no spatial spread)
r0 = 0;
pearsonR0 = 0;
r023 = 0;
pearsonR023 = 0;
eqPrev = 0;
eqPrevStd = 0;
r0Adap = 0;
pearsonR0Adap = 0;
r0Adap23 = 0;
pearsonR0Adap23 = 0;
eqPrevAdap = 0;
eqPrevStdAdap = 0;
r0Sigmoid = 0;
r0Sigmoid23 = 0;
eqPrevSigmoid = 0;

waveSpeed = 0;
pearsonWaveSpeed = 0;
info = 0;

fixTR0 = 2;
fixTMin = 3;
fixTMax = 9;

r0Factor = 4;
averBeg = 5;
averEnd = 12;

% Simulation
[tPrev, prevR0, prevR023, popR0, prevEP, popEP, kMonitored, tDistance, maxDistance, proportionReached] = ...
    Main(xInit, yInit, infectionRate, 1/latentPeriod, contactRate, diagnoseProp, diagnosePeriod, migrationRate, K);

% Check if disease could spread
if ( max(maxDistance) < minDistance && proportionReached < minProportionReached)
    % Check if it could spread at least locally
    if ( max(prevR0) > minPrev )
        proportionReached = -1;
        return;
    else
        return;
    end
end

% Get wave speed
if (info ~= 1)
    [waveSpeed, pearsonWaveSpeed] = GetWaveSpeed(tDistance, maxDistance);
end

% Get r0, eqPrev, duration
%    Fixed period
%      r0
[r0, pearsonR0] = GetR0Fix(tPrev, prevR0, fixTR0, kMonitored);
[r023, pearsonR023] = GetR0Fix(tPrev, prevR023, fixTR0, kMonitored);
% eqPrev
[eqPrev, eqPrevStd] = GetEqPrevFixed(tPrev, prevEP, popEP, fixTMin, fixTMax);


%    Adaptive period - r0
[r0Adap, pearsonR0Adap] = GetR0Adap(tPrev, prevR0, popR0, kMonitored, r0Factor);
[r0Adap23, pearsonR0Adap23] = GetR0Adap(tPrev, prevR023, popR0, kMonitored, r0Factor);
%   Adaptive period - eqPrev
[eqPrevAdap, eqPrevStdAdap] = GetEqPrevAdap(tPrev, prevEP, popEP, kMonitored, r0Factor, averBeg, averEnd);

% r0Sigmoid, eqPrevSigmoid
% From 1 diseased animal to 30 animal, no zeros
[r0Sigmoid, temp1, temp2] = GetSigmoid(tPrev, prevR0, popR0);
[r0Sigmoid23, temp1, temp2] = GetSigmoid(tPrev, prevR023, popR0);
[temp1, temp2, eqPrevSigmoid] = GetSigmoid(tPrev, prevEP, popEP);
