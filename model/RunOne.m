% bestSims = importdata('../../FitModel/ABC/best200_BestSimsParamStats_Obs0.txt','\t',1);
% 
% i=127;
% xInit = bestSims.data(i,9);
% yInit = bestSims.data(i,10);
% infectionConst = bestSims.data(i,3);
% outbreakConst = 1/bestSims.data(i,4);
% globalTransmissionConst = bestSims.data(i,5);
% diagnoseProp = bestSims.data(i,7);
% diagnosePeriod = bestSims.data(i,8);
% migrationConst = bestSims.data(i,6);
% popAll=60000;
% rng(1);

% [extinctionTime, relPop, meanPrev, timeRet, pop, prev, nInf, popInit, prevInit, occupied, ...
%     popSpatial, prevSpatial, nInfSpatial, popAge, nInfAge] = MainLong(xInit, ...
%     yInit, infectionConst, outbreakConst, globalTransmissionConst, diagnoseProp, ...
%     diagnosePeriod, migrationConst, popAll);

infectionConst = 30.2090;
outbreakConst = 1/0.2262;
globalTransmissionConst = 0.8129;
migrationConst = 0.3929;
diagnoseProp = 0.9319;
diagnosePeriod = 0.4706;
xInit = 9; 
yInit = 23;
popAll = 60000;

[tPrev, prevR0, prevR023, popR0, popR023, prevEP, popEP, tDistance, maxDistance, proportionReached] = ...
    MainFit(xInit, yInit, infectionConst, outbreakConst, globalTransmissionConst, ...
    diagnoseProp, diagnosePeriod, migrationConst, popAll)