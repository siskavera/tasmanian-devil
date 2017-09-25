function [prop, population] = GetEqDist(deathRate, birthRate, ageTransferConst)

nAgeClasses = length(deathRate);
if ( nargin < 3 )
    ageTransferConst = ones(1,nAgeClasses);
end

population = zeros(size(deathRate));
alpha = population;

alpha(1) = 1;
for i = 2:(nAgeClasses - 1)
    alpha(i) = alpha(i-1)*ageTransferConst(i-1)/(ageTransferConst(i) + deathRate(i));
end
alpha(end) = alpha(end-1)*ageTransferConst(end-1)/deathRate(end);

population(1) = (1 - (ageTransferConst(1) + deathRate(1))/sum(birthRate.*alpha))/sum(alpha);
population = alpha.*population(1);

prop = sum(population);
population = population/sum(population);