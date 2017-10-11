function RunLongTerms(nSims, runID, iSimBeg)
%  matlab -nodisplay -nojvm -r "RunLongTerms(200,1,1);exit;"
rng('shuffle');
popAll = 60000; % From McCallum, personal communication

if nargin < 3
    iSimBeg = 1;
end

bestSims = importdata('best200_BestSimsParamStats_Obs0.txt','\t',1);
nSims = min([size( bestSims.data, 1), nSims]);
[~, ordering] = sort( bestSims.data(:,2), 'ascend' );
bestSims.data = bestSims.data(ordering, :);

for i = iSimBeg:nSims
    xInit = bestSims.data(i,9);
    yInit = bestSims.data(i,10);
    infectionConst = bestSims.data(i,3);
    outbreakConst = 1/bestSims.data(i,4);
    globalTransmissionConst = bestSims.data(i,5);
    diagnoseProp = bestSims.data(i,7);
    diagnosePeriod = bestSims.data(i,8);
    migrationConst = bestSims.data(i,6);
    
    fprintf('Running simulation no. %i\n', i);
    fprintf('  Origin:                 [%i,%i]\n', xInit, yInit);
    fprintf('  Infection rate:         %.4f\n', infectionConst);
    fprintf('  Latent period           %.4f\n', 1/outbreakConst);
    fprintf('  Diagnosable period:     %.4f\n', diagnosePeriod);
    fprintf('  Diagnosable proportion: %.4f\n', diagnoseProp);
    fprintf('  Contact rate:           %.4f\n', globalTransmissionConst);
    fprintf('  Migration rate:         %.4f\n\n', migrationConst);
    
    [extinctionTime, relPop, meanPrev, timeRet, pop, prev, nInf, popInit, prevInit, occupied, ...
        popSpatial, prevSpatial, nInfSpatial, popAge, nInfAge] = ...
        MainLong(xInit, yInit, infectionConst, outbreakConst, globalTransmissionConst, ...
        diagnoseProp, diagnosePeriod, migrationConst, popAll);
    
    fprintf('  Extinction time:        %.4f\n', extinctionTime);
    fprintf('  Relative population:    %.4f\n', relPop);
    fprintf('  Mean prevalence:        %.4f\n\n\n', meanPrev);
    
    save(sprintf('../scratch/LongTerm%d/result_%d', runID, i), '-v7.3');
end
