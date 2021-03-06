function Sampler(nSimuMax, outputID)
% SAMPLER Evaluate Monte Carlo samples
%   Sampler(nSimuMax, outputID) runs nSimuMax simulations and saves
%   parameters and summary statistics to result_{outputID}.mat
%
% Parameter ranges hard-coded
%   param   min max
%   eta     0   50
%   c       0   1
%   m       0   2
%   prop    0   1
%   diaPer  0   1.2
%   latPer  0   1
%   yInit   1   42
%   xInit   1   42
%
% Saved variables are:
%
%   para Table of input parameters
%     eta    Local infection rate
%     latPer Latent period
%     c      Between-patch contact rate
%     m      Migration rate
%     prop   Proportion of diagnosable period when animals are not yet
%     infectious
%     diaPer Total diagnosable period
%     xInit  x coordinate of origin (first mutation)
%     yInit  y coordinate of origin (first mutation)
%
%   data Table of summary statistics
%
% r0, pearsonR0, r023, pearsonR023, eqPrev, eqPrevStd, ...
%        r0Adap, pearsonR0Adap, r0Adap23, pearsonR0Adap23, eqPrevAdap, eqPrevStdAdap, ...
%        r0Sigmoid, r0Sigmoid23, eqPrevSigmoid, waveSpeed, pearsonWaveSpeed, ...
%        waveSpeedTheo, waveSpeedTheoPearson, proportionReached

% Seed random number generator
rng('shuffle')

popAll = 60000; % From McCallum, personal communication
K = ReadK(100000); % Doesn't matter, it will only be used to determine if an origin is valid
isSea = (K < 20);

[waveSpeedTheoAll, waveSpeedTheoPearsonAll] = ReadWaveSpeedTheo();

% All parameters are drawn from a uniform distribution in [parMin:parMax]
etaMin = 0;         % Infection rate
etaMax = 50;
cMin = 0;           % Damping factor for betwepropen-patch contact rate
cMax = 1;
mMin = 0;           % Migration rate
mMax = 2;
propMin = 0;        % tD/(tD+tI)
propMax = 1;
diaPerMin = 0;   % tD+tI
diaPerMax = 1.2;
latPerMin = 0;   % tE
latPerMax = 1;
yInitMin = 1;
yInitMax = 42;
xInitMin = 1;
xInitMax = 42;

data = zeros(nSimuMax,20);
para = zeros(nSimuMax,8);

N = 100;

for i = 1:nSimuMax
    eta    = etaMin    + (etaMax    -    etaMin) * rand;
    c      = cMin      + (cMax      -      cMin) * rand;
    m      = mMin      + (mMax      -      mMin) * rand;
    prop   = propMin   + (propMax   -   propMin) * rand;
    diaPer = diaPerMin + (diaPerMax - diaPerMin) * rand;
    latPer = latPerMin + (latPerMax - latPerMin) * rand;
    
    xInit  = floor(xInitMin  + (xInitMax  -  xInitMin + 1) * rand);
    yInit  = floor(yInitMin  + (yInitMax  -  yInitMin + 1) * rand);
    while ( isSea(xInit,yInit) )
        xInit  = floor(xInitMin  + (xInitMax  -  xInitMin + 1) * rand);
        yInit  = floor(yInitMin  + (yInitMax  -  yInitMin + 1) * rand);
    end
    
    [r0, pearsonR0, r023, pearsonR023, eqPrev, eqPrevStd, ...
        r0Adap, pearsonR0Adap, r0Adap23, pearsonR0Adap23, eqPrevAdap, eqPrevStdAdap, ...
        r0Sigmoid, r0Sigmoid23, eqPrevSigmoid, waveSpeed, pearsonWaveSpeed, proportionReached] = ...
        GetSummaryStats(eta, latPer, c, prop, diaPer, m, xInit, yInit, popAll);
     waveSpeedTheo = waveSpeedTheoAll(xInit, yInit);
     waveSpeedTheoPearson = waveSpeedTheoPearsonAll(xInit, yInit);
     
    fprintf('%d %.6f %.6f %.6f %.6f %.6f %.6f %d %d \n', ...
        i, eta, c, m, prop, diaPer, m, xInit, yInit);
    fprintf('%d %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n\n', ...
        i, r0, eqPrev,  r0Adap, eqPrevAdap, r0Sigmoid, eqPrevSigmoid, waveSpeed);
    
    data(i,:) = [r0, pearsonR0, r023, pearsonR023, eqPrev, eqPrevStd, ...
        r0Adap, pearsonR0Adap, r0Adap23, pearsonR0Adap23, eqPrevAdap, eqPrevStdAdap, ...
        r0Sigmoid, r0Sigmoid23, eqPrevSigmoid, waveSpeed, pearsonWaveSpeed, ...
        waveSpeedTheo, waveSpeedTheoPearson, proportionReached];
    para(i,:) = [eta, latPer, c, m, prop, diaPer, xInit, yInit];
    
    if (mod(i,5) == 0)
        save(sprintf('result_%d', outputID), 'para', 'data');
    end
end

save(sprintf('result_%d', outputID), 'para', 'data');
