function [eqPrev, eqPrevStd] = GetEqPrevFixed(t, prev, pop, tMin, tMax)
% GETEQPREVFIXED Calculate equilibrial prevalence from a fixed period
%   [eqPrev, eqPrevStd] = GetEqPrevFixed(t, prev, pop, tMin, tMax) returns 
%   the the mean (eqPrev) and standard deviation (eqPrevStd) of equilibrial
%   prevalence.

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

%      eqPrev
temp = [];
for i = 1:nEqPrev
    thisPrev = prevEPAll{i};
    thisPrev = thisPrev(1:10:end);   % Considering only quarterly field trips
    thisT = tPrevEP{i};
    thisT = thisT(1:10:end);
    temp = [temp; thisPrev(thisT > tMin & thisT < tMax)];
end

if (size(temp,1) < 2)
    eqPrev = 0;
    eqPrevStd = 0;
else
    eqPrev = mean(temp);
    eqPrevStd = std(temp);
end