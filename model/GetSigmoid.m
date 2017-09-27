function [r0, tHalf, eqPrev] = GetSigmoid(t, prev, pop)
% GETSIGMOID Calculate summary statistics from a sigmoid fit
%   [r0, tHalf, eqPrev] = GetSigmoid(t, prev, pop) returns fitted initial
%   increase rate r0, half-life tHalf and equilibrial value eqPrev from a
%   sigmoidal fit to prevalence (prev) as a function of time (t). Only
%   considers time-series until population size drops below 30, to avoid
%   highly stochastic period

nSites = size(prev,2);
tPrevEP = cell(1,nSites);
prevEPAll = cell(1,nSites);
popEPAll = cell(1,nSites);

for i = 1:nSites
    iMin = find( prev(:,i)>0, 1);
    tPrevEP{i} = t(iMin:end)-t(iMin);
    prevEPAll{i} = prev(iMin:end,i);
    popEPAll{i} = pop(iMin:end,i);
end


if ( ~isempty(find(pop<30,1 )) > 0 )
    iEnd = find(pop < 30, 1);
else
    iEnd = length(pop);
end

if ( isempty(iEnd) || iEnd<20 )
    eqPrev = 0;
    tHalf = 0;
    r0 = 0;
    return;
end

eqPrev = 0;
r0 = 0;
tHalf = 0;
nValid = 0;
for i = 1:nSites
    if ( ~isempty(find(popEPAll{i}<30,1 )) > 0 )
        iEnd = find(popEPAll{i} < 30, 1);
    else
        iEnd = length(popEPAll{i});
    end
    if ( ~isempty(iEnd) && iEnd>=20 )
        try
            [temp1, thistHalf, thisEqPrev, thisR0] = fit_logistic( tPrevEP{i}(1:10:iEnd), prevEPAll{i}(1:10:iEnd) );
        catch
            thistHalf = nan;
            thisEqPrev = nan;
            thisR0 = nan;
        end
        if ( ~isnan(thistHalf) && ~isnan(thisEqPrev) && ~isnan(thisR0))
            tHalf = tHalf + thistHalf;
            eqPrev = eqPrev + thisEqPrev;
            r0 = r0 + thisR0;
            nValid = nValid + 1;
        end
    end
end

if nValid == 0
    eqPrev = 0;
    tHalf = 0;
    r0 = 0;
end
