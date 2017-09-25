function [x,y] = ChoosePatch(rate)

n = size(rate,1);
vectorVersion = rate(:);
cumsumVector = cumsum(vectorVersion)/sum(vectorVersion);
index = find(cumsumVector > rand,1);
y = ceil(index/n);
x = index - (y-1)*n;