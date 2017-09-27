function [r0, pearsonR] = GetR0Fix(t, prev, fixT, kMonitored)
% GETR0FIX Calculate initial increase in prevalence from a fixed period
%   [r0, pearsonR] = GetR0Fix(t, prev, fixT, kMonitored) returns estimated
%   initial increase rate in prevalence (r0) and the pearson correlation
%   coefficient from the fit (pearsonR)
%
% r0 is estimated from a linear fit to the logarithm of prevalence. The
% first fixT years after the first case and quarterly field trips are 
% considered. 0.1/kMonitored is added to the prevalence to avoid problems 
% from disease disappearing while fitting. This does not influence the 
% result of the fit much.

iMin = find( prev>0, 1);
t = t(iMin:end) - t(iMin);
prev = prev(iMin:end);

%    Fixed period - Consider every 10th point only!
%      r0

temp = [t(t < fixT), prev(t < fixT)+0.1/kMonitored];
temp = temp(1:10:end,:);
if (size(temp,1) < 3)
    r0 = 0;
    pearsonR = 0;
    info = 3;
else
    p = polyfit(temp(:,1), log(temp(:,2)), 1);
    r0 = p(1);
    pearsonR = corrcoef(temp(:,1), log(temp(:,2)));
    pearsonR = pearsonR(1,end);
end

