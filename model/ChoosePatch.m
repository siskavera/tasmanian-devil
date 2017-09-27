function [x,y] = ChoosePatch(rate)
% CHOOSEPATCH Choose patch for a an event
%   ChoosePatch(rate) Returns indices for the chosen patch of the event,
%   chosen with probability proportional to its rate, given the rate matrix
n = size(rate,1);
vectorVersion = rate(:);
cumsumVector = cumsum(vectorVersion)/sum(vectorVersion);
index = find(cumsumVector > rand,1);
y = ceil(index/n);
x = index - (y-1)*n;