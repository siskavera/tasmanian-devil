function [r0, pearsonR] = GetR0Adap(t, prev, pop, prevForPop, r0Factor)
% GETR0ADAP Calculate initial increase in prevalence from an adaptive period
%   [r0, pearsonR] = GetR0Adap(t, prev, pop, r0Factor)
%
% Adaptive period: r0 calculated from 3 diseased animals until prevalence
% surpasses maximal prevalence / r0Factor

iBeg = find(prevForPop.*pop >= 3, 1);
iEnd = find(prev>(max(prev)/r0Factor),1);

if ( isempty(iEnd) || isempty(iBeg) || iEnd<(iBeg+19))
    r0 = 0;
    pearsonR = 0;
else
    temp = [t(iBeg:iEnd), prev(iBeg:iEnd)];
    temp = temp(1:10:end,:);
    p = polyfit(temp(:,1), log(temp(:,2)), 1);
    r0 = p(1);
    pearsonR = corrcoef(temp(:,1), log(temp(:,2)));
    pearsonR = pearsonR(1,end);
end
