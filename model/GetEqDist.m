function [prop, population] = GetEqDist(deathRate, birthRate, ageTransferConst)
% GETEQDIST Calculates equilibrial age distribution
%   [prop, population] = GetEqDist(deathRate, birthRate,
%   ageTransferConst)returns equilibrial relative population size (prop)
%   and age distribution (population), given deathRate, birthRate and
%   transfer rate between age classes.

nAgeClasses = length(deathRate);
if ( nargin < 3 )
    ageTransferConst = ones(1,nAgeClasses); % Yearly age-classes
end

population = zeros(size(deathRate));
alpha = population;

alpha(1) = 1;
for i = 2:(nAgeClasses - 1)
    % Ratio between age classes converges to ratio of input and output
    % rates
    alpha(i) = alpha(i-1)*ageTransferConst(i-1)/(ageTransferConst(i) + deathRate(i));
end
alpha(end) = alpha(end-1)*ageTransferConst(end-1)/deathRate(end);

% From input = output for p0
%   birth rate     carrying capacity  death   aging
% sum_i(p_i b_i) * (1 - sum_i(p_i)) - (d_0 + age_0) * p_0 = 0
% and p_i = alpha_i * p_0
population(1) = (1 - (ageTransferConst(1) + deathRate(1))/sum(birthRate.*alpha) ) / sum(alpha);
population = alpha.*population(1);

prop = sum(population);
population = population/sum(population);