function [eqPrev, eqPrevStd] = GetEqPrevAdap(t, prev, pop, kMonitored, r0Factor, averBeg, averEnd)
% GETEQPREVADAP Calculate equilibrial prevalence from an adaptive period
%   [eqPrev, eqPrevStd] = GetEqPrevAdap(t, prev, pop, kMonitored, r0Factor,
%   averBeg, averEnd) returns the mean (eqPrev) and standard deviation
%   (eqPrevStd) of equilibrial prevalence.
%
% Adaptive period: r0 calculated from 5 diseased animals until prevalence
% surpasses maximal prevalence / r0Factor. Equibrial prevalence then
% calculated from time > averBeg/r0 until time > averEnd/r0

nEqPrev = size(prev,2);
tPrevEP = cell(1,nEqPrev);
prevEPAll = cell(1,nEqPrev);
popEPAll = cell(1,nEqPrev);

for i = 1:nEqPrev
    iMin = find( prev(:,i)>0, 1);
    tPrevEP{i} = t(iMin:end)-t(iMin);
    prevEPAll{i} = prev(iMin:end,i);
    popEPAll{i} = pop(iMin:end,i);
end

r0EP = 0;
nGood = 0;
for i = 1:nEqPrev
    thisT = tPrevEP{i};
    thisPrev = prevEPAll{i};
    thisPop = popEPAll{i};

    iBeg = find(thisPrev.*thisPop > 5, 1);
    iEnd = find(thisPrev>(max(thisPrev)/r0Factor),1);
    if ~( isempty(iEnd) || isempty(iBeg) || iEnd<(iBeg+19))
        temp = [thisT(iBeg:iEnd), thisPrev(iBeg:iEnd)+0.1/kMonitored];
        temp = temp(1:10:end,:);
        p = polyfit(temp(:,1), log(temp(:,2)), 1);
        r0EP = r0EP + p(1);
        nGood = nGood + 1;
    end
end

r0EP = r0EP/nGood;
prevEPLate = [];

for i = 1:nEqPrev
    thisT = tPrevEP{i};
    thisPrev = prevEPAll{i};
    
    iBeg = find(thisT > averBeg/r0EP, 1);
    iEnd = find(thisT > averEnd/r0EP, 1);
    if ( ~isempty(iEnd) && ~isempty(iBeg) )
        prevEPLate = [prevEPLate; thisPrev(iBeg:10:iEnd)];
    end
end

if ( length(prevEPLate) > 2 )
    eqPrev = mean(prevEPLate);
    eqPrevStd = std(prevEPLate);
    
else
    eqPrev = 0;
    eqPrevStd = 0;
end