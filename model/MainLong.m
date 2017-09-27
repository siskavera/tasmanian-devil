function [extinctionTime, relPop, meanPrev, timeRet, pop, prev, nInf, popInit, prevInit, occupied, ...
    popSpatial, prevSpatial, nInfSpatial, popAge, nInfAge] = ...
    MainLong(xInit, yInit, infectionConst, outbreakConst, globalTransmissionConst, ...
    diagnoseProp, diagnosePeriod, migrationConst, K)
% Continuous time, stochastic (finite, discrete number of individuals)
% spatial model with Poisson birth and death processes, patches coupled by
% migration and nearest neighbour contact, age structure

nSteps = 1000000000;
tTransient = 4;
tMax = 250;
nInit = 2;

% Setting up spatial environment
population = K;
[n, m] = size(population);
isSea = (K==0);
K(K == 0) = inf;

% Indices of initially infected patch
% R0: Fentonbury [24,22] ~75
% eqPrev: Mt William [7,34] 200-217, Bronte Park [19,21] ~70, Buckland
% [24,30] ??
monitoredR0 = [24 22];
monitoredEqPrev = [7, 34; 19, 21; 24, 30];
nMonitoredEqPrev = size(monitoredEqPrev,1);
kMonitored = K(monitoredR0(1), monitoredR0(2) );
kMonitoredEP = zeros(1, nMonitoredEqPrev);
for i = 1:nMonitoredEqPrev
    kMonitoredEP = K(monitoredEqPrev(i,1), monitoredEqPrev(i,2));
end
prevInit = nInit/K(xInit, yInit);

% Setting parameters
nAgeClasses = 5;
birthConst = [0 0.13 1.3 1.65 1.21]';
deathConst = [0.2526 0.2439 0.1214 0.0719 0.8146]'; % From Beeton, scaled
% deathConst = [0.211 0.144 0.133 0.146 0.353];  % From McCallum
diagnoseConst = 1/(diagnoseProp*diagnosePeriod);
removalConst = 1/((1-diagnoseProp)*diagnosePeriod);

migConst = [0 migrationConst migrationConst migrationConst migrationConst]'; % Juveniles assumed not to migrate
transmissionConst = [0 0           0           0           0; ...
    0 0.602*infectionConst 0.602*infectionConst 0.602*infectionConst 0.602*infectionConst; ...
    0 0.602*infectionConst infectionConst       infectionConst       infectionConst; ...
    0 0.602*infectionConst infectionConst       infectionConst       infectionConst; ...
    0 0.602*infectionConst infectionConst       infectionConst       infectionConst]';

[eqPopProp, eqAgeDist] = GetEqDist(deathConst, birthConst);  % From deterministic dynamics
temp = population;
for i = 1:nAgeClasses
    population(:,:,i) = ceil(temp*eqAgeDist(i)*eqPopProp);
end
kMonitored = kMonitored * eqPopProp;
kMonitoredEP = kMonitoredEP * eqPopProp;

xSize = 410;
ySize = 410;

% Initial state - healthy
susceptible = population;
exposed = zeros(n,m,nAgeClasses);
diagnosable = zeros(n,m,nAgeClasses);
infective = zeros(n,m,nAgeClasses);

% Assume a band of uninhabited patches on the borders
infectiveNeighbour = zeros(n,m,nAgeClasses);

% Constructing rates
birthRate = zeros(n,m);
DMRateS = zeros(n,m);  % Death, aging and migration together
DMRateE = zeros(n,m);
DMRateD = zeros(n,m);
DMRateI = zeros(n,m);
transmissionRate = zeros(n,m);

for i = 1:n
    for j = 1:m
        birthRate(i,j) = birthConst'*reshape(population(i,j,:),nAgeClasses,1) * max( [0, 1-sum(population(i,j,:))/K(i,j)]);
        DMRateS(i,j) = sum((deathConst + migConst + 1) .* reshape(susceptible(i,j,:),nAgeClasses,1));
        DMRateE(i,j) = sum((deathConst + migConst + outbreakConst + 1) .* reshape(exposed(i,j,:),nAgeClasses,1));
        DMRateD(i,j) = sum((deathConst + migConst + diagnoseConst + 1) .* reshape(diagnosable(i,j,:),nAgeClasses,1));
        DMRateI(i,j) = sum((deathConst + migConst + removalConst + 1)  .* reshape(infective(i,j,:),nAgeClasses,1));
        transmissionRate(i,j) = (reshape(infective(i,j,:),1,nAgeClasses) + globalTransmissionConst*reshape(infectiveNeighbour(i,j,:),1,nAgeClasses))...
            *transmissionConst* reshape(susceptible(i,j,:),nAgeClasses,1)/sum(population(i,j,:));
        transmissionRate(isnan(transmissionRate)) = 0;
    end
end

time = 0;
rates = [sum(sum(DMRateS,2),1) sum(sum(DMRateE,2),1) sum(sum(DMRateD,2),1) sum(sum(DMRateI,2),1)...
        sum(sum(transmissionRate,2),1) sum(sum(birthRate,1),2)];
    
for iStep = 1:nSteps
    
    rateSum = sum(rates);
    if (rateSum == 0)
        break;
    end
    
    time = time - log(rand)/rateSum;
    
    r1 = rand;
    r2 = rand;
    r3 = rand;
    eventNo = find(cumsum(rates)/rateSum > r1, 1);
    
    % Indicator of migration
    isMig = 0;
    
    switch eventNo
        case 1 % Susceptible: death, migration or aging
            [x,y] = ChoosePatch(DMRateS);
            sOnPatch = sum(susceptible(x,y,:));
            
            if (r3 < sOnPatch/DMRateS(x,y)) % Aging
                ageClass = find(cumsum(susceptible(x,y,:))/sum(susceptible(x,y,:)) > r2, 1);
                newAge = min([nAgeClasses, ageClass+1]);
                susceptible(x,y,newAge) = susceptible(x,y,newAge) + 1;
                population(x,y,newAge) = population(x,y,newAge) + 1;
                susceptible(x,y,ageClass) = susceptible(x,y,ageClass) - 1;
                population(x,y,ageClass) = population(x,y,ageClass) - 1;
                
            else if (r3 < sOnPatch*(1+migrationConst*sum(susceptible(x,y,2:end))/sOnPatch)/DMRateS(x,y)) % Migration
                    ageClass = find(cumsum(susceptible(x,y,2:end))/sum(susceptible(x,y,2:end)) > r2, 1)+1;
                    [newX,newY] = ChoosePatchForMigrant(isSea, x, y);
                    susceptible(newX,newY,ageClass) = susceptible(newX,newY,ageClass) + 1;
                    population(newX,newY,ageClass) = population(newX,newY,ageClass) + 1;
                    susceptible(x,y,ageClass) = susceptible(x,y,ageClass) - 1;
                    population(x,y,ageClass) = population(x,y,ageClass) - 1;
                    
                    thisPop = sum(population(newX,newY,:));
                    if (thisPop > 0)
                        new = (reshape(infective(newX,newY,:),1,nAgeClasses) + ...
                            globalTransmissionConst*reshape(infectiveNeighbour(newX,newY,:),1,nAgeClasses))...
                            *transmissionConst* reshape(susceptible(newX,newY,:),nAgeClasses,1)/thisPop;
                        rates(5) = rates(5) - transmissionRate(newX,newY) + new;
                        transmissionRate(newX,newY) = new;
                    else
                        rates(5) = rates(5) - transmissionRate(newX, newY);
                        transmissionRate(newX,newY) = 0;
                    end
                    new = birthConst'*...
                        reshape(population(newX,newY,:),nAgeClasses,1) * max( [0, 1-sum(population(newX,newY,:))/K(newX,newY)]);
                    rates(6) = rates(6) - birthRate(newX,newY) + new;
                    birthRate(newX,newY) = new;
                    new = sum((deathConst + migConst + 1) .* reshape(susceptible(newX,newY,:),nAgeClasses,1));
                    rates(1) = rates(1) - DMRateS(newX,newY) + new;
                    DMRateS(newX,newY) = new;
                else % Death
                    ageClass = find(cumsum(deathConst.*reshape(susceptible(x,y,:),nAgeClasses,1))/sum(deathConst.*reshape(susceptible(x,y,:),nAgeClasses,1)) > r2, 1);
                    susceptible(x,y,ageClass) = susceptible(x,y,ageClass) - 1;
                    population(x,y,ageClass) = population(x,y,ageClass) - 1;
                end
            end
            
            new = birthConst'*reshape(population(x,y,:),nAgeClasses,1) * max( [0, 1-sum(population(x,y,:))/K(x,y)]);
            rates(6) = rates(6) - birthRate(x,y) + new;
            birthRate(x,y) = new;
            new = sum((deathConst + migConst + 1) .* reshape(susceptible(x,y,:),nAgeClasses,1));
            rates(1) = rates(1) - DMRateS(x,y) + new;
            DMRateS(x,y) = new;
            thisPop = sum(population(x,y,:));
            if (thisPop > 0)
                new = (reshape(infective(x,y,:),1,nAgeClasses) + ...
                    globalTransmissionConst*reshape(infectiveNeighbour(x,y,:),1,nAgeClasses))...
                    *transmissionConst* reshape(susceptible(x,y,:),nAgeClasses,1)/thisPop;
                rates(5) = rates(5) - transmissionRate(x,y) + new;
                transmissionRate(x,y) = new;
            else
                rates(5) = rates(5) - transmissionRate(x,y);
                transmissionRate(x,y) = 0;
            end
            
        case 2 % Exposed: death, migration, aging or outbreak
            [x,y] = ChoosePatch(DMRateE);
            eOnPatch = sum(exposed(x,y,:));
            
            if (r3 < eOnPatch/DMRateE(x,y)) % Aging
                ageClass = find(cumsum(exposed(x,y,:))/sum(exposed(x,y,:)) > r2, 1);
                newAge = min([nAgeClasses, ageClass+1]);
                exposed(x,y,newAge) = exposed(x,y,newAge) + 1;
                population(x,y,newAge) = population(x,y,newAge) + 1;
                exposed(x,y,ageClass) = exposed(x,y,ageClass) - 1;
                population(x,y,ageClass) = population(x,y,ageClass) - 1;
                
            else if (r3 < eOnPatch*(1+migrationConst*sum(exposed(x,y,2:end))/eOnPatch)/DMRateE(x,y)) % Migration
                    ageClass = find(cumsum(exposed(x,y,2:end))/sum(exposed(x,y,2:end)) > r2, 1)+1;
                    [newX,newY] = ChoosePatchForMigrant(isSea, x, y);
                    exposed(newX,newY,ageClass) = exposed(newX,newY,ageClass) + 1;
                    population(newX,newY,ageClass) = population(newX,newY,ageClass) + 1;
                    exposed(x,y,ageClass) = exposed(x,y,ageClass) - 1;
                    population(x,y,ageClass) = population(x,y,ageClass) - 1;
                    
                    thisPop = sum(population(newX,newY,:));
                    if (thisPop > 0)
                        new = (reshape(infective(newX,newY,:),1,nAgeClasses) + ...
                            globalTransmissionConst*reshape(infectiveNeighbour(newX,newY,:),1,nAgeClasses))...
                            *transmissionConst* reshape(susceptible(newX,newY,:),nAgeClasses,1)/thisPop;
                        rates(5) = rates(5) - transmissionRate(newX,newY) + new;
                        transmissionRate(newX,newY) = new;
                    else
                        rates(5) = rates(5) - transmissionRate(newX,newY);
                        transmissionRate(newX,newY) = 0;
                    end
                    
                    new = birthConst'*...
                        reshape(population(newX,newY,:),nAgeClasses,1) * max( [0, 1-sum(population(newX,newY,:))/K(newX,newY)]);
                    rates(6) = rates(6) - birthRate(newX,newY) + new;
                    birthRate(newX,newY) = new;
                    new = sum((deathConst + migConst + outbreakConst +  1) .* reshape(exposed(newX,newY,:),nAgeClasses,1));
                    rates(2) = rates(2) - DMRateE(newX,newY) + new;
                    DMRateE(newX,newY) = new;
                    
                else if (r3 < eOnPatch*(1+migrationConst*sum(exposed(x,y,2:end))/eOnPatch+outbreakConst)/DMRateE(x,y)) % Outbreak
                        ageClass = find(cumsum(exposed(x,y,:))/sum(exposed(x,y,:)) > r2, 1);
                        diagnosable(x,y,ageClass) = diagnosable(x,y,ageClass) + 1;
                        exposed(x,y,ageClass) = exposed(x,y,ageClass) - 1;
                        new = sum((deathConst + migConst + diagnoseConst + 1) .* reshape(diagnosable(x,y,:),nAgeClasses,1));
                        rates(3) = rates(3) - DMRateD(x,y) + new;
                        DMRateD(x,y) = new;
                    else % Death
                        ageClass = find(cumsum(deathConst.*reshape(exposed(x,y,:),nAgeClasses,1))/sum(deathConst.*reshape(exposed(x,y,:),nAgeClasses,1)) > r2, 1);
                        exposed(x,y,ageClass) = exposed(x,y,ageClass) - 1;
                        population(x,y,ageClass) = population(x,y,ageClass) - 1;
                    end
                end
            end
            
            new = birthConst'*reshape(population(x,y,:),nAgeClasses,1) * max( [0, 1-sum(population(x,y,:))/K(x,y)]);
            rates(6) = rates(6) - birthRate(x,y) + new;
            birthRate(x,y) = new;
            new = sum((deathConst + migConst + outbreakConst + 1) .* reshape(exposed(x,y,:),nAgeClasses,1));
            rates(2) = rates(2) - DMRateE(x,y) + new;
            DMRateE(x,y) = new;
            thisPop = sum(population(x,y,:));
            if (thisPop > 0)
                new = (reshape(infective(x,y,:),1,nAgeClasses) + ...
                    globalTransmissionConst*reshape(infectiveNeighbour(x,y,:),1,nAgeClasses))...
                    *transmissionConst* reshape(susceptible(x,y,:),nAgeClasses,1)/thisPop;
                rates(5) = rates(5) - transmissionRate(x,y) + new;
                transmissionRate(x,y) = new;
            else
                rates(5) = rates(5) - transmissionRate(x,y);
                transmissionRate(x,y) = 0;
            end
            
        case 3 % Diagnosable: death, migration, aging or outbreak
            [x,y] = ChoosePatch(DMRateD);
            dOnPatch = sum(diagnosable(x,y,:));
            
            if (r3 < dOnPatch/DMRateD(x,y)) % Aging
                ageClass = find(cumsum(diagnosable(x,y,:))/sum(diagnosable(x,y,:)) > r2, 1);
                newAge = min([nAgeClasses, ageClass+1]);
                diagnosable(x,y,newAge) = diagnosable(x,y,newAge) + 1;
                population(x,y,newAge) = population(x,y,newAge) + 1;
                diagnosable(x,y,ageClass) = diagnosable(x,y,ageClass) - 1;
                population(x,y,ageClass) = population(x,y,ageClass) - 1;
                
            else if (r3 < dOnPatch*(1+migrationConst*sum(diagnosable(x,y,2:end))/dOnPatch)/DMRateD(x,y)) % Migration
                    ageClass = find(cumsum(diagnosable(x,y,2:end))/sum(diagnosable(x,y,2:end)) > r2, 1) + 1;
                    [newX,newY] = ChoosePatchForMigrant(isSea, x, y);
                    diagnosable(newX,newY,ageClass) = diagnosable(newX,newY,ageClass) + 1;
                    population(newX,newY,ageClass) = population(newX,newY,ageClass) + 1;
                    diagnosable(x,y,ageClass) = diagnosable(x,y,ageClass) - 1;
                    population(x,y,ageClass) = population(x,y,ageClass) - 1;
                    
                    thisPop = sum(population(newX,newY,:));
                    if (thisPop > 0)
                        new = (reshape(infective(newX,newY,:),1,nAgeClasses) + ...
                            globalTransmissionConst*reshape(infectiveNeighbour(newX,newY,:),1,nAgeClasses))...
                            *transmissionConst* reshape(susceptible(newX,newY,:),nAgeClasses,1)/thisPop;
                        rates(5) = rates(5) - transmissionRate(newX,newY) + new;
                        transmissionRate(newX,newY) = new;
                    else
                        rates(5) = rates(5) - transmissionRate(newX,newY);
                        transmissionRate(newX,newY) = 0;
                    end
                    
                   new = birthConst'*...
                       reshape(population(newX,newY,:),nAgeClasses,1) * max( [0, 1-sum(population(newX,newY,:))/K(newX,newY)]);
                   rates(6) = rates(6) -  birthRate(newX,newY) + new;
                   birthRate(newX,newY)= new;
                   new = sum((deathConst + migConst + diagnoseConst +  1) .* reshape(diagnosable(newX,newY,:),nAgeClasses,1));
                   rates(3) = rates(3) - DMRateD(newX,newY) + new;
                   DMRateD(newX,newY) = new;
                   
                else if (r3 < dOnPatch*(1+migrationConst*sum(diagnosable(x,y,2:end))/dOnPatch+diagnoseConst)/DMRateD(x,y)) % Diagnose
                        ageClass = find(cumsum(diagnosable(x,y,:))/sum(diagnosable(x,y,:)) > r2, 1);
                        infective(x,y,ageClass) = infective(x,y,ageClass) + 1;
                        diagnosable(x,y,ageClass) = diagnosable(x,y,ageClass) - 1;
                        infectiveNeighbour(x+1,y,ageClass) = infectiveNeighbour(x+1,y,ageClass) + 1;
                        infectiveNeighbour(x-1,y,ageClass) = infectiveNeighbour(x-1,y,ageClass) + 1;
                        infectiveNeighbour(x,y+1,ageClass) = infectiveNeighbour(x,y+1,ageClass) + 1;
                        infectiveNeighbour(x,y-1,ageClass) = infectiveNeighbour(x,y-1,ageClass) + 1;
                        
                        indices = [ [x-1; x-1; x-1; x; x; x; x+1; x+1; x+1] repmat([y-1; y; y+1],3,1)];
                        for i = 1:size(indices,1);
                            thisX = indices(i,1);
                            thisY = indices(i,2);
                            thisPop = sum(population(thisX,thisY,:));
                            if (thisPop > 0)
                                new = (reshape(infective(thisX,thisY,:),1,nAgeClasses) + ...
                                    globalTransmissionConst*reshape(infectiveNeighbour(thisX,thisY,:),1,nAgeClasses))...
                                    *transmissionConst* reshape(susceptible(thisX,thisY,:),nAgeClasses,1)/thisPop;
                                rates(5) = rates(5) - transmissionRate(thisX,thisY) + new;
                                transmissionRate(thisX,thisY) = new;
                            else
                                rates(5) = rates(5) - transmissionRate(thisX,thisY);
                                transmissionRate(thisX,thisY) = 0;
                            end
                        end
                        
                        new = sum((deathConst + migConst + removalConst + 1) .* reshape(infective(x,y,:),nAgeClasses,1));
                        rates(4) = rates(4) - DMRateI(x,y) + new;
                        DMRateI(x,y) = new;
                    else % Death
                        ageClass = find(cumsum(deathConst.*reshape(diagnosable(x,y,:),nAgeClasses,1))/sum(deathConst.*reshape(diagnosable(x,y,:),nAgeClasses,1)) > r2, 1);
                        diagnosable(x,y,ageClass) = diagnosable(x,y,ageClass) - 1;
                        population(x,y,ageClass) = population(x,y,ageClass) - 1;
                    end
                end
            end
            
            new = birthConst'*reshape(population(x,y,:),nAgeClasses,1) * max( [0, 1-sum(population(x,y,:))/K(x,y)]);
            rates(6) = rates(6) - birthRate(x,y) + new;
            birthRate(x,y) = new;
            new = sum((deathConst + migConst + diagnoseConst + 1) .* reshape(diagnosable(x,y,:),nAgeClasses,1));
            rates(3) = rates(3) - DMRateD(x,y) + new;
            DMRateD(x,y) = new;
            thisPop = sum(population(x,y,:));
            if (thisPop > 0)
                new = (reshape(infective(x,y,:),1,nAgeClasses) + ...
                    globalTransmissionConst*reshape(infectiveNeighbour(x,y,:),1,nAgeClasses))...
                    *transmissionConst* reshape(susceptible(x,y,:),nAgeClasses,1)/thisPop;
                rates(5) = rates(5) - transmissionRate(x,y) + new;
                transmissionRate(x,y) = new;
            else
                rates(5) = rates(5) - transmissionRate(x,y);
                transmissionRate(x,y) = 0;
            end
            
        case 4 % Infective: death, migration or aging
            [x,y] = ChoosePatch(DMRateI);
            iOnPatch = sum(infective(x,y,:));
            
            if (r3 < iOnPatch/DMRateI(x,y)) % Aging
                ageClass = find(cumsum(infective(x,y,:))/sum(infective(x,y,:)) > r2, 1);
                newAge = min([nAgeClasses, ageClass+1]);
                infective(x,y,newAge) = infective(x,y,newAge) + 1;
                population(x,y,newAge) = population(x,y,newAge) + 1;
                infectiveNeighbour(x+1,y,ageClass) = infectiveNeighbour(x+1,y,ageClass) - 1;
                infectiveNeighbour(x-1,y,ageClass) = infectiveNeighbour(x-1,y,ageClass) - 1;
                infectiveNeighbour(x,y+1,ageClass) = infectiveNeighbour(x,y+1,ageClass) - 1;
                infectiveNeighbour(x,y-1,ageClass) = infectiveNeighbour(x,y-1,ageClass) - 1;
                infectiveNeighbour(x+1,y,newAge) = infectiveNeighbour(x+1,y,newAge) + 1;
                infectiveNeighbour(x-1,y,newAge) = infectiveNeighbour(x-1,y,newAge) + 1;
                infectiveNeighbour(x,y+1,newAge) = infectiveNeighbour(x,y+1,newAge) + 1;
                infectiveNeighbour(x,y-1,newAge) = infectiveNeighbour(x,y-1,newAge) + 1;
                infective(x,y,ageClass) = infective(x,y,ageClass) - 1;
                population(x,y,ageClass) = population(x,y,ageClass) - 1;
                
            else if (r3 < iOnPatch*(1+migrationConst*sum(infective(x,y,2:end))/iOnPatch)/DMRateI(x,y)) % Migration
                    ageClass = find(cumsum(infective(x,y,2:end))/sum(infective(x,y,2:end)) > r2, 1)+1;
                    [newX,newY] = ChoosePatchForMigrant(isSea, x, y);
                    infective(newX,newY,ageClass) = infective(newX,newY,ageClass) + 1;
                    population(newX,newY,ageClass) = population(newX,newY,ageClass) + 1;
                    infectiveNeighbour(x+1,y,ageClass) = infectiveNeighbour(x+1,y,ageClass) - 1;
                    infectiveNeighbour(x-1,y,ageClass) = infectiveNeighbour(x-1,y,ageClass) - 1;
                    infectiveNeighbour(x,y+1,ageClass) = infectiveNeighbour(x,y+1,ageClass) - 1;
                    infectiveNeighbour(x,y-1,ageClass) = infectiveNeighbour(x,y-1,ageClass) - 1;
                    infectiveNeighbour(newX+1,newY,ageClass) = infectiveNeighbour(newX+1,newY,ageClass) + 1;
                    infectiveNeighbour(newX-1,newY,ageClass) = infectiveNeighbour(newX-1,newY,ageClass) + 1;
                    infectiveNeighbour(newX,newY+1,ageClass) = infectiveNeighbour(newX,newY+1,ageClass) + 1;
                    infectiveNeighbour(newX,newY-1,ageClass) = infectiveNeighbour(newX,newY-1,ageClass) + 1;
                    infective(x,y,ageClass) = infective(x,y,ageClass) - 1;
                    population(x,y,ageClass) = population(x,y,ageClass) - 1;
                    
                    isMig = 1;
                    
                    new = birthConst'*...
                        reshape(population(newX,newY,:),nAgeClasses,1) * max( [0, 1-sum(population(newX,newY,:))/K(newX,newY)]);
                    rates(6) = rates(6) - birthRate(newX,newY) + new;
                    birthRate(newX,newY) = new;
                    new = sum((deathConst + migConst + removalConst + 1) .* reshape(infective(newX,newY,:),nAgeClasses,1));
                    rates(4) = rates(4) - DMRateI(newX,newY) + new;
                    DMRateI(newX,newY) = new;
                    
                else % Death
                    ageClass = find(cumsum(reshape(infective(x,y,:),nAgeClasses,1).*(deathConst+removalConst))/...
                        sum(reshape(infective(x,y,:),nAgeClasses,1).*(deathConst+removalConst)) > r2, 1);
                    infectiveNeighbour(x+1,y,ageClass) = infectiveNeighbour(x+1,y,ageClass) - 1;
                    infectiveNeighbour(x-1,y,ageClass) = infectiveNeighbour(x-1,y,ageClass) - 1;
                    infectiveNeighbour(x,y+1,ageClass) = infectiveNeighbour(x,y+1,ageClass) - 1;
                    infectiveNeighbour(x,y-1,ageClass) = infectiveNeighbour(x,y-1,ageClass) - 1;
                    infective(x,y,ageClass) = infective(x,y,ageClass) - 1;
                    population(x,y,ageClass) = population(x,y,ageClass) - 1;
                end
            end
            
            indices = [ [x-1; x-1; x-1; x; x; x; x+1; x+1; x+1] repmat([y-1; y; y+1],3,1)];
            if (isMig == 1)
                newIndices = [ [newX-1; newX-1; newX-1; newX; newX; newX; newX+1; newX+1; newX+1] ...
                    repmat([newY-1; newY; newY+1],3,1)];
                indices = [indices; newIndices];
            end
            
            for i = 1:size(indices,1);
                thisX = indices(i,1);
                thisY = indices(i,2);
                thisPop = sum(population(thisX,thisY,:));
                if (thisPop > 0)
                    new = (reshape(infective(thisX,thisY,:),1,nAgeClasses) + ...
                        globalTransmissionConst*reshape(infectiveNeighbour(thisX,thisY,:),1,nAgeClasses))...
                        *transmissionConst* reshape(susceptible(thisX,thisY,:),nAgeClasses,1)/thisPop;
                    rates(5) = rates(5) - transmissionRate(thisX,thisY) + new;
                    transmissionRate(thisX, thisY) = new;
                else
                    rates(5) = rates(5) - transmissionRate(thisX,thisY);
                    transmissionRate(thisX,thisY) = 0;
                end
            end
            
            new = birthConst'*reshape(population(x,y,:),nAgeClasses,1) * max( [0, 1-sum(population(x,y,:))/K(x,y)]);
            rates(6) = rates(6) - birthRate(x,y) + new;
            birthRate(x,y) = new;
            new = sum((deathConst + migConst + removalConst + 1) .* reshape(infective(x,y,:),nAgeClasses,1));
            rates(4) = rates(4) - DMRateI(x,y) + new;
            DMRateI(x,y) = new;
            
            thisPop = sum(population(x,y,:));
            if (thisPop > 0)
                new = (reshape(infective(x,y,:),1,nAgeClasses) + ...
                    globalTransmissionConst*reshape(infectiveNeighbour(x,y,:),1,nAgeClasses))...
                    *transmissionConst* reshape(susceptible(x,y,:),nAgeClasses,1)/thisPop;
                rates(5) = rates(5) - transmissionRate(x,y) + new;
                transmissionRate(x,y) = new;
            else
                rates(5) = rates(5) - transmissionRate(x,y);
                transmissionRate(x,y) = 0;
            end
            
        case 5 % Transmission
            [x,y] = ChoosePatch(transmissionRate);
            temp = ( (reshape(infective(x,y,:),1,nAgeClasses) + ...
                globalTransmissionConst*reshape(infectiveNeighbour(x,y,:),1,nAgeClasses))...
                * transmissionConst ) .* reshape(susceptible(x,y,:),1,nAgeClasses);
            ageClass = find(cumsum(temp)/sum(temp) > r2, 1);
            susceptible(x,y,ageClass) = susceptible(x,y,ageClass) - 1;
            exposed(x,y,ageClass) = exposed(x,y,ageClass) + 1;
            
            new = sum((deathConst + migConst + 1) .* reshape(susceptible(x,y,:),nAgeClasses,1));
            rates(1) = rates(1) - DMRateS(x,y) + new;
            DMRateS(x,y) = new;
            new = sum((deathConst + migConst + outbreakConst + 1) .* reshape(exposed(x,y,:),nAgeClasses,1));
            rates(2) = rates(2) - DMRateE(x,y) + new;
            DMRateE(x,y) = new;
            
            thisPop = sum(population(x,y,:));
            if (thisPop > 0)
                new = (reshape(infective(x,y,:),1,nAgeClasses) + ...
                    globalTransmissionConst*reshape(infectiveNeighbour(x,y,:),1,nAgeClasses))...
                    *transmissionConst* reshape(susceptible(x,y,:),nAgeClasses,1)/thisPop;
                rates(5) = rates(5) - transmissionRate(x,y) + new;
                transmissionRate(x,y)= new;
            else
                rates(5) = rates(5) - transmissionRate(x,y);
                transmissionRate(x,y) = 0;
            end

        case 6 % Birth
            [x,y] = ChoosePatch(birthRate);
            population(x,y,1) = population(x,y,1) + 1;
            susceptible(x,y,1) = susceptible(x,y,1) + 1;
            
            new = birthConst'*reshape(population(x,y,:),nAgeClasses,1) * max( [0, 1-sum(population(x,y,:))/K(x,y)]);
            rates(6) = rates(6) - birthRate(x,y) + new;
            birthRate(x,y) = new;
            new = sum((deathConst + migConst + 1) .* reshape(susceptible(x,y,:),nAgeClasses,1));
            rates(1) = rates(1) - DMRateS(x,y) + new;
            DMRateS(x,y) = new;
            
            thisPop = sum(population(x,y,:));
            if (thisPop > 0)
                new = (reshape(infective(x,y,:),1,nAgeClasses) + ...
                    globalTransmissionConst*reshape(infectiveNeighbour(x,y,:),1,nAgeClasses))...
                    *transmissionConst* reshape(susceptible(x,y,:),nAgeClasses,1)/thisPop;
                rates(5) = rates(5) - transmissionRate(x,y) + new;
                transmissionRate(x,y) = new;
                
            else
                rates(5) = rates(5) - transmissionRate(x,y);
                transmissionRate(x,y) = 0;
            end      
    end
    

    if (time > tTransient)       
        break;
    end
end
fprintf('%d\n', sum(sum(sum(population))));
time = 0;

% Infection
infective(xInit, yInit, 2:end) = ceil(prevInit*population(xInit, yInit, 2:end));
susceptible(xInit, yInit, 2:end) = susceptible(xInit, yInit, 2:end) - ceil(prevInit*population(xInit, yInit, 2:end));

% Setting up monitored variables
isInfected = zeros(n,m);

timeRet = [];
pop = [];
prev = [];
popInit = [];
prevInit = [];
occupied = [];
popSpatial = [];
prevSpatial = [];
nInfSpatial = [];
nInf = [];
popAge = [];
nInfAge = [];

iRecord = 1;
tMonitored = 0;

extinctionTime = -1;
relPop = 0;
meanPrev = 0;

isInfected(xInit, yInit) = 1;

infectiveNeighbour(2:end-1,2:end-1,:) = infective(1:end-2,2:end-1,:) + infective(3:end,2:end-1,:) +...
    infective(2:end-1,1:end-2,:) + infective(2:end-1,3:end,:);

% Constructing rates
birthRate = zeros(n,m);
DMRateS = zeros(n,m);  % Death, aging and migration together
DMRateE = zeros(n,m);
DMRateD = zeros(n,m);
DMRateI = zeros(n,m);
transmissionRate = zeros(n,m);

for i = 1:n
    for j = 1:m
        birthRate(i,j) = birthConst'*reshape(population(i,j,:),nAgeClasses,1) * max( [0, 1-sum(population(i,j,:))/K(i,j)]);
        DMRateS(i,j) = sum((deathConst + migConst + 1) .* reshape(susceptible(i,j,:),nAgeClasses,1));
        DMRateE(i,j) = sum((deathConst + migConst + outbreakConst + 1) .* reshape(exposed(i,j,:),nAgeClasses,1));
        DMRateD(i,j) = sum((deathConst + migConst + diagnoseConst + 1) .* reshape(diagnosable(i,j,:),nAgeClasses,1));
        DMRateI(i,j) = sum((deathConst + migConst + removalConst + 1)  .* reshape(infective(i,j,:),nAgeClasses,1));
        transmissionRate(i,j) = (reshape(infective(i,j,:),1,nAgeClasses) + globalTransmissionConst*reshape(infectiveNeighbour(i,j,:),1,nAgeClasses))...
            *transmissionConst* reshape(susceptible(i,j,:),nAgeClasses,1)/sum(population(i,j,:));
        transmissionRate(isnan(transmissionRate)) = 0;
    end
end

rates = [sum(sum(DMRateS,2),1) sum(sum(DMRateE,2),1) sum(sum(DMRateD,2),1) sum(sum(DMRateI,2),1)...
        sum(sum(transmissionRate,2),1) sum(sum(birthRate,1),2)];
for iStep = 1:nSteps
    
    rateSum = sum(rates);
    if (rateSum == 0)
        break;
    end
    
    time = time - log(rand)/rateSum;
    
    r1 = rand;
    r2 = rand;
    r3 = rand;
    eventNo = find(cumsum(rates)/rateSum > r1, 1);
    
    % Indicator of migration
    isMig = 0;
    
    switch eventNo
        case 1 % Susceptible: death, migration or aging
            [x,y] = ChoosePatch(DMRateS);
            sOnPatch = sum(susceptible(x,y,:));
            
            if (r3 < sOnPatch/DMRateS(x,y)) % Aging
                ageClass = find(cumsum(susceptible(x,y,:))/sum(susceptible(x,y,:)) > r2, 1);
                newAge = min([nAgeClasses, ageClass+1]);
                susceptible(x,y,newAge) = susceptible(x,y,newAge) + 1;
                population(x,y,newAge) = population(x,y,newAge) + 1;
                susceptible(x,y,ageClass) = susceptible(x,y,ageClass) - 1;
                population(x,y,ageClass) = population(x,y,ageClass) - 1;
                
            else if (r3 < sOnPatch*(1+migrationConst*sum(susceptible(x,y,2:end))/sOnPatch)/DMRateS(x,y)) % Migration
                    ageClass = find(cumsum(susceptible(x,y,2:end))/sum(susceptible(x,y,2:end)) > r2, 1)+1;
                    [newX,newY] = ChoosePatchForMigrant(isSea, x, y);
                    susceptible(newX,newY,ageClass) = susceptible(newX,newY,ageClass) + 1;
                    population(newX,newY,ageClass) = population(newX,newY,ageClass) + 1;
                    susceptible(x,y,ageClass) = susceptible(x,y,ageClass) - 1;
                    population(x,y,ageClass) = population(x,y,ageClass) - 1;
                    
                    thisPop = sum(population(newX,newY,:));
                    if (thisPop > 0)
                        new = (reshape(infective(newX,newY,:),1,nAgeClasses) + ...
                            globalTransmissionConst*reshape(infectiveNeighbour(newX,newY,:),1,nAgeClasses))...
                            *transmissionConst* reshape(susceptible(newX,newY,:),nAgeClasses,1)/thisPop;
                        rates(5) = rates(5) - transmissionRate(newX,newY) + new;
                        transmissionRate(newX,newY) = new;
                    else
                        rates(5) = rates(5) - transmissionRate(newX, newY);
                        transmissionRate(newX,newY) = 0;
                    end
                    new = birthConst'*...
                        reshape(population(newX,newY,:),nAgeClasses,1) * max( [0, 1-sum(population(newX,newY,:))/K(newX,newY)]);
                    rates(6) = rates(6) - birthRate(newX,newY) + new;
                    birthRate(newX,newY) = new;
                    new = sum((deathConst + migConst + 1) .* reshape(susceptible(newX,newY,:),nAgeClasses,1));
                    rates(1) = rates(1) - DMRateS(newX,newY) + new;
                    DMRateS(newX,newY) = new;
                else % Death
                    ageClass = find(cumsum(deathConst.*reshape(susceptible(x,y,:),nAgeClasses,1))/sum(deathConst.*reshape(susceptible(x,y,:),nAgeClasses,1)) > r2, 1);
                    susceptible(x,y,ageClass) = susceptible(x,y,ageClass) - 1;
                    population(x,y,ageClass) = population(x,y,ageClass) - 1;
                end
            end
            
            new = birthConst'*reshape(population(x,y,:),nAgeClasses,1) * max( [0, 1-sum(population(x,y,:))/K(x,y)]);
            rates(6) = rates(6) - birthRate(x,y) + new;
            birthRate(x,y) = new;
            new = sum((deathConst + migConst + 1) .* reshape(susceptible(x,y,:),nAgeClasses,1));
            rates(1) = rates(1) - DMRateS(x,y) + new;
            DMRateS(x,y) = new;
            thisPop = sum(population(x,y,:));
            if (thisPop > 0)
                new = (reshape(infective(x,y,:),1,nAgeClasses) + ...
                    globalTransmissionConst*reshape(infectiveNeighbour(x,y,:),1,nAgeClasses))...
                    *transmissionConst* reshape(susceptible(x,y,:),nAgeClasses,1)/thisPop;
                rates(5) = rates(5) - transmissionRate(x,y) + new;
                transmissionRate(x,y) = new;
            else
                rates(5) = rates(5) - transmissionRate(x,y);
                transmissionRate(x,y) = 0;
            end
            
        case 2 % Exposed: death, migration, aging or outbreak
            [x,y] = ChoosePatch(DMRateE);
            eOnPatch = sum(exposed(x,y,:));
            
            if (r3 < eOnPatch/DMRateE(x,y)) % Aging
                ageClass = find(cumsum(exposed(x,y,:))/sum(exposed(x,y,:)) > r2, 1);
                newAge = min([nAgeClasses, ageClass+1]);
                exposed(x,y,newAge) = exposed(x,y,newAge) + 1;
                population(x,y,newAge) = population(x,y,newAge) + 1;
                exposed(x,y,ageClass) = exposed(x,y,ageClass) - 1;
                population(x,y,ageClass) = population(x,y,ageClass) - 1;
                
            else if (r3 < eOnPatch*(1+migrationConst*sum(exposed(x,y,2:end))/eOnPatch)/DMRateE(x,y)) % Migration
                    ageClass = find(cumsum(exposed(x,y,2:end))/sum(exposed(x,y,2:end)) > r2, 1)+1;
                    [newX,newY] = ChoosePatchForMigrant(isSea, x, y);
                    exposed(newX,newY,ageClass) = exposed(newX,newY,ageClass) + 1;
                    population(newX,newY,ageClass) = population(newX,newY,ageClass) + 1;
                    exposed(x,y,ageClass) = exposed(x,y,ageClass) - 1;
                    population(x,y,ageClass) = population(x,y,ageClass) - 1;
                    
                    thisPop = sum(population(newX,newY,:));
                    if (thisPop > 0)
                        new = (reshape(infective(newX,newY,:),1,nAgeClasses) + ...
                            globalTransmissionConst*reshape(infectiveNeighbour(newX,newY,:),1,nAgeClasses))...
                            *transmissionConst* reshape(susceptible(newX,newY,:),nAgeClasses,1)/thisPop;
                        rates(5) = rates(5) - transmissionRate(newX,newY) + new;
                        transmissionRate(newX,newY) = new;
                    else
                        rates(5) = rates(5) - transmissionRate(newX,newY);
                        transmissionRate(newX,newY) = 0;
                    end
                    
                    new = birthConst'*...
                        reshape(population(newX,newY,:),nAgeClasses,1) * max( [0, 1-sum(population(newX,newY,:))/K(newX,newY)]);
                    rates(6) = rates(6) - birthRate(newX,newY) + new;
                    birthRate(newX,newY) = new;
                    new = sum((deathConst + migConst + outbreakConst +  1) .* reshape(exposed(newX,newY,:),nAgeClasses,1));
                    rates(2) = rates(2) - DMRateE(newX,newY) + new;
                    DMRateE(newX,newY) = new;
                    
                else if (r3 < eOnPatch*(1+migrationConst*sum(exposed(x,y,2:end))/eOnPatch+outbreakConst)/DMRateE(x,y)) % Outbreak
                        ageClass = find(cumsum(exposed(x,y,:))/sum(exposed(x,y,:)) > r2, 1);
                        diagnosable(x,y,ageClass) = diagnosable(x,y,ageClass) + 1;
                        exposed(x,y,ageClass) = exposed(x,y,ageClass) - 1;
                        new = sum((deathConst + migConst + diagnoseConst + 1) .* reshape(diagnosable(x,y,:),nAgeClasses,1));
                        rates(3) = rates(3) - DMRateD(x,y) + new;
                        DMRateD(x,y) = new;
                    else % Death
                        ageClass = find(cumsum(deathConst.*reshape(exposed(x,y,:),nAgeClasses,1))/sum(deathConst.*reshape(exposed(x,y,:),nAgeClasses,1)) > r2, 1);
                        exposed(x,y,ageClass) = exposed(x,y,ageClass) - 1;
                        population(x,y,ageClass) = population(x,y,ageClass) - 1;
                    end
                end
            end
            
            new = birthConst'*reshape(population(x,y,:),nAgeClasses,1) * max( [0, 1-sum(population(x,y,:))/K(x,y)]);
            rates(6) = rates(6) - birthRate(x,y) + new;
            birthRate(x,y) = new;
            new = sum((deathConst + migConst + outbreakConst + 1) .* reshape(exposed(x,y,:),nAgeClasses,1));
            rates(2) = rates(2) - DMRateE(x,y) + new;
            DMRateE(x,y) = new;
            thisPop = sum(population(x,y,:));
            if (thisPop > 0)
                new = (reshape(infective(x,y,:),1,nAgeClasses) + ...
                    globalTransmissionConst*reshape(infectiveNeighbour(x,y,:),1,nAgeClasses))...
                    *transmissionConst* reshape(susceptible(x,y,:),nAgeClasses,1)/thisPop;
                rates(5) = rates(5) - transmissionRate(x,y) + new;
                transmissionRate(x,y) = new;
            else
                rates(5) = rates(5) - transmissionRate(x,y);
                transmissionRate(x,y) = 0;
            end
            
        case 3 % Diagnosable: death, migration, aging or outbreak
            [x,y] = ChoosePatch(DMRateD);
            dOnPatch = sum(diagnosable(x,y,:));
            
            if (r3 < dOnPatch/DMRateD(x,y)) % Aging
                ageClass = find(cumsum(diagnosable(x,y,:))/sum(diagnosable(x,y,:)) > r2, 1);
                newAge = min([nAgeClasses, ageClass+1]);
                diagnosable(x,y,newAge) = diagnosable(x,y,newAge) + 1;
                population(x,y,newAge) = population(x,y,newAge) + 1;
                diagnosable(x,y,ageClass) = diagnosable(x,y,ageClass) - 1;
                population(x,y,ageClass) = population(x,y,ageClass) - 1;
                
            else if (r3 < dOnPatch*(1+migrationConst*sum(diagnosable(x,y,2:end))/dOnPatch)/DMRateD(x,y)) % Migration
                    ageClass = find(cumsum(diagnosable(x,y,2:end))/sum(diagnosable(x,y,2:end)) > r2, 1) + 1;
                    [newX,newY] = ChoosePatchForMigrant(isSea, x, y);
                    diagnosable(newX,newY,ageClass) = diagnosable(newX,newY,ageClass) + 1;
                    population(newX,newY,ageClass) = population(newX,newY,ageClass) + 1;
                    diagnosable(x,y,ageClass) = diagnosable(x,y,ageClass) - 1;
                    population(x,y,ageClass) = population(x,y,ageClass) - 1;
                    
                    thisPop = sum(population(newX,newY,:));
                    if (thisPop > 0)
                        new = (reshape(infective(newX,newY,:),1,nAgeClasses) + ...
                            globalTransmissionConst*reshape(infectiveNeighbour(newX,newY,:),1,nAgeClasses))...
                            *transmissionConst* reshape(susceptible(newX,newY,:),nAgeClasses,1)/thisPop;
                        rates(5) = rates(5) - transmissionRate(newX,newY) + new;
                        transmissionRate(newX,newY) = new;
                    else
                        rates(5) = rates(5) - transmissionRate(newX,newY);
                        transmissionRate(newX,newY) = 0;
                    end
                    
                   new = birthConst'*...
                       reshape(population(newX,newY,:),nAgeClasses,1) * max( [0, 1-sum(population(newX,newY,:))/K(newX,newY)]);
                   rates(6) = rates(6) -  birthRate(newX,newY) + new;
                   birthRate(newX,newY)= new;
                   new = sum((deathConst + migConst + diagnoseConst +  1) .* reshape(diagnosable(newX,newY,:),nAgeClasses,1));
                   rates(3) = rates(3) - DMRateD(newX,newY) + new;
                   DMRateD(newX,newY) = new;
                   
                else if (r3 < dOnPatch*(1+migrationConst*sum(diagnosable(x,y,2:end))/dOnPatch+diagnoseConst)/DMRateD(x,y)) % Diagnose
                        ageClass = find(cumsum(diagnosable(x,y,:))/sum(diagnosable(x,y,:)) > r2, 1);
                        infective(x,y,ageClass) = infective(x,y,ageClass) + 1;
                        diagnosable(x,y,ageClass) = diagnosable(x,y,ageClass) - 1;
                        infectiveNeighbour(x+1,y,ageClass) = infectiveNeighbour(x+1,y,ageClass) + 1;
                        infectiveNeighbour(x-1,y,ageClass) = infectiveNeighbour(x-1,y,ageClass) + 1;
                        infectiveNeighbour(x,y+1,ageClass) = infectiveNeighbour(x,y+1,ageClass) + 1;
                        infectiveNeighbour(x,y-1,ageClass) = infectiveNeighbour(x,y-1,ageClass) + 1;
                        
                        indices = [ [x-1; x-1; x-1; x; x; x; x+1; x+1; x+1] repmat([y-1; y; y+1],3,1)];
                        for i = 1:size(indices,1);
                            thisX = indices(i,1);
                            thisY = indices(i,2);
                            thisPop = sum(population(thisX,thisY,:));
                            if (thisPop > 0)
                                new = (reshape(infective(thisX,thisY,:),1,nAgeClasses) + ...
                                    globalTransmissionConst*reshape(infectiveNeighbour(thisX,thisY,:),1,nAgeClasses))...
                                    *transmissionConst* reshape(susceptible(thisX,thisY,:),nAgeClasses,1)/thisPop;
                                rates(5) = rates(5) - transmissionRate(thisX,thisY) + new;
                                transmissionRate(thisX,thisY) = new;
                            else
                                rates(5) = rates(5) - transmissionRate(thisX,thisY);
                                transmissionRate(thisX,thisY) = 0;
                            end
                        end
                        
                        new = sum((deathConst + migConst + removalConst + 1) .* reshape(infective(x,y,:),nAgeClasses,1));
                        rates(4) = rates(4) - DMRateI(x,y) + new;
                        DMRateI(x,y) = new;
                    else % Death
                        ageClass = find(cumsum(deathConst.*reshape(diagnosable(x,y,:),nAgeClasses,1))/sum(deathConst.*reshape(diagnosable(x,y,:),nAgeClasses,1)) > r2, 1);
                        diagnosable(x,y,ageClass) = diagnosable(x,y,ageClass) - 1;
                        population(x,y,ageClass) = population(x,y,ageClass) - 1;
                    end
                end
            end
            
            new = birthConst'*reshape(population(x,y,:),nAgeClasses,1) * max( [0, 1-sum(population(x,y,:))/K(x,y)]);
            rates(6) = rates(6) - birthRate(x,y) + new;
            birthRate(x,y) = new;
            new = sum((deathConst + migConst + diagnoseConst + 1) .* reshape(diagnosable(x,y,:),nAgeClasses,1));
            rates(3) = rates(3) - DMRateD(x,y) + new;
            DMRateD(x,y) = new;
            thisPop = sum(population(x,y,:));
            if (thisPop > 0)
                new = (reshape(infective(x,y,:),1,nAgeClasses) + ...
                    globalTransmissionConst*reshape(infectiveNeighbour(x,y,:),1,nAgeClasses))...
                    *transmissionConst* reshape(susceptible(x,y,:),nAgeClasses,1)/thisPop;
                rates(5) = rates(5) - transmissionRate(x,y) + new;
                transmissionRate(x,y) = new;
            else
                rates(5) = rates(5) - transmissionRate(x,y);
                transmissionRate(x,y) = 0;
            end
            
        case 4 % Infective: death, migration or aging
            [x,y] = ChoosePatch(DMRateI);
            iOnPatch = sum(infective(x,y,:));
            
            if (r3 < iOnPatch/DMRateI(x,y)) % Aging
                ageClass = find(cumsum(infective(x,y,:))/sum(infective(x,y,:)) > r2, 1);
                newAge = min([nAgeClasses, ageClass+1]);
                infective(x,y,newAge) = infective(x,y,newAge) + 1;
                population(x,y,newAge) = population(x,y,newAge) + 1;
                infectiveNeighbour(x+1,y,ageClass) = infectiveNeighbour(x+1,y,ageClass) - 1;
                infectiveNeighbour(x-1,y,ageClass) = infectiveNeighbour(x-1,y,ageClass) - 1;
                infectiveNeighbour(x,y+1,ageClass) = infectiveNeighbour(x,y+1,ageClass) - 1;
                infectiveNeighbour(x,y-1,ageClass) = infectiveNeighbour(x,y-1,ageClass) - 1;
                infectiveNeighbour(x+1,y,newAge) = infectiveNeighbour(x+1,y,newAge) + 1;
                infectiveNeighbour(x-1,y,newAge) = infectiveNeighbour(x-1,y,newAge) + 1;
                infectiveNeighbour(x,y+1,newAge) = infectiveNeighbour(x,y+1,newAge) + 1;
                infectiveNeighbour(x,y-1,newAge) = infectiveNeighbour(x,y-1,newAge) + 1;
                infective(x,y,ageClass) = infective(x,y,ageClass) - 1;
                population(x,y,ageClass) = population(x,y,ageClass) - 1;
                
            else if (r3 < iOnPatch*(1+migrationConst*sum(infective(x,y,2:end))/iOnPatch)/DMRateI(x,y)) % Migration
                    ageClass = find(cumsum(infective(x,y,2:end))/sum(infective(x,y,2:end)) > r2, 1)+1;
                    [newX,newY] = ChoosePatchForMigrant(isSea, x, y);
                    infective(newX,newY,ageClass) = infective(newX,newY,ageClass) + 1;
                    population(newX,newY,ageClass) = population(newX,newY,ageClass) + 1;
                    infectiveNeighbour(x+1,y,ageClass) = infectiveNeighbour(x+1,y,ageClass) - 1;
                    infectiveNeighbour(x-1,y,ageClass) = infectiveNeighbour(x-1,y,ageClass) - 1;
                    infectiveNeighbour(x,y+1,ageClass) = infectiveNeighbour(x,y+1,ageClass) - 1;
                    infectiveNeighbour(x,y-1,ageClass) = infectiveNeighbour(x,y-1,ageClass) - 1;
                    infectiveNeighbour(newX+1,newY,ageClass) = infectiveNeighbour(newX+1,newY,ageClass) + 1;
                    infectiveNeighbour(newX-1,newY,ageClass) = infectiveNeighbour(newX-1,newY,ageClass) + 1;
                    infectiveNeighbour(newX,newY+1,ageClass) = infectiveNeighbour(newX,newY+1,ageClass) + 1;
                    infectiveNeighbour(newX,newY-1,ageClass) = infectiveNeighbour(newX,newY-1,ageClass) + 1;
                    infective(x,y,ageClass) = infective(x,y,ageClass) - 1;
                    population(x,y,ageClass) = population(x,y,ageClass) - 1;
                    
                    isMig = 1;
                    
                    new = birthConst'*...
                        reshape(population(newX,newY,:),nAgeClasses,1) * max( [0, 1-sum(population(newX,newY,:))/K(newX,newY)]);
                    rates(6) = rates(6) - birthRate(newX,newY) + new;
                    birthRate(newX,newY) = new;
                    new = sum((deathConst + migConst + removalConst + 1) .* reshape(infective(newX,newY,:),nAgeClasses,1));
                    rates(4) = rates(4) - DMRateI(newX,newY) + new;
                    DMRateI(newX,newY) = new;
                    
                else % Death
                    ageClass = find(cumsum(reshape(infective(x,y,:),nAgeClasses,1).*(deathConst+removalConst))/...
                        sum(reshape(infective(x,y,:),nAgeClasses,1).*(deathConst+removalConst)) > r2, 1);
                    infectiveNeighbour(x+1,y,ageClass) = infectiveNeighbour(x+1,y,ageClass) - 1;
                    infectiveNeighbour(x-1,y,ageClass) = infectiveNeighbour(x-1,y,ageClass) - 1;
                    infectiveNeighbour(x,y+1,ageClass) = infectiveNeighbour(x,y+1,ageClass) - 1;
                    infectiveNeighbour(x,y-1,ageClass) = infectiveNeighbour(x,y-1,ageClass) - 1;
                    infective(x,y,ageClass) = infective(x,y,ageClass) - 1;
                    population(x,y,ageClass) = population(x,y,ageClass) - 1;
                end
            end
            
            indices = [ [x-1; x-1; x-1; x; x; x; x+1; x+1; x+1] repmat([y-1; y; y+1],3,1)];
            if (isMig == 1)
                newIndices = [ [newX-1; newX-1; newX-1; newX; newX; newX; newX+1; newX+1; newX+1] ...
                    repmat([newY-1; newY; newY+1],3,1)];
                indices = [indices; newIndices];
            end
            
            for i = 1:size(indices,1);
                thisX = indices(i,1);
                thisY = indices(i,2);
                thisPop = sum(population(thisX,thisY,:));
                if (thisPop > 0)
                    new = (reshape(infective(thisX,thisY,:),1,nAgeClasses) + ...
                        globalTransmissionConst*reshape(infectiveNeighbour(thisX,thisY,:),1,nAgeClasses))...
                        *transmissionConst* reshape(susceptible(thisX,thisY,:),nAgeClasses,1)/thisPop;
                    rates(5) = rates(5) - transmissionRate(thisX,thisY) + new;
                    transmissionRate(thisX,thisY) = new;
                else
                    rates(5) = rates(5) - transmissionRate(thisX,thisY);
                    transmissionRate(thisX,thisY) = 0;
                end
            end
            
            new = birthConst'*reshape(population(x,y,:),nAgeClasses,1) * max( [0, 1-sum(population(x,y,:))/K(x,y)]);
            rates(6) = rates(6) - birthRate(x,y) + new;
            birthRate(x,y) = new;
            new = sum((deathConst + migConst + removalConst + 1) .* reshape(infective(x,y,:),nAgeClasses,1));
            rates(4) = rates(4) - DMRateI(x,y) + new;
            DMRateI(x,y) = new;
            
            thisPop = sum(population(x,y,:));
            if (thisPop > 0)
                new = (reshape(infective(x,y,:),1,nAgeClasses) + ...
                    globalTransmissionConst*reshape(infectiveNeighbour(x,y,:),1,nAgeClasses))...
                    *transmissionConst* reshape(susceptible(x,y,:),nAgeClasses,1)/thisPop;
                rates(5) = rates(5) - transmissionRate(x,y) + new;
                transmissionRate(x,y) = new;
            else
                rates(5) = rates(5) - transmissionRate(x,y);
                transmissionRate(x,y) = 0;
            end
            
        case 5 % Transmission
            [x,y] = ChoosePatch(transmissionRate);
            temp = ( (reshape(infective(x,y,:),1,nAgeClasses) + ...
                globalTransmissionConst*reshape(infectiveNeighbour(x,y,:),1,nAgeClasses))...
                * transmissionConst ) .* reshape(susceptible(x,y,:),1,nAgeClasses);
            ageClass = find(cumsum(temp)/sum(temp) > r2, 1);
            susceptible(x,y,ageClass) = susceptible(x,y,ageClass) - 1;
            exposed(x,y,ageClass) = exposed(x,y,ageClass) + 1;
            
            new = sum((deathConst + migConst + 1) .* reshape(susceptible(x,y,:),nAgeClasses,1));
            rates(1) = rates(1) - DMRateS(x,y) + new;
            DMRateS(x,y) = new;
            new = sum((deathConst + migConst + outbreakConst + 1) .* reshape(exposed(x,y,:),nAgeClasses,1));
            rates(2) = rates(2) - DMRateE(x,y) + new;
            DMRateE(x,y) = new;
            
            thisPop = sum(population(x,y,:));
            if (thisPop > 0)
                new = (reshape(infective(x,y,:),1,nAgeClasses) + ...
                    globalTransmissionConst*reshape(infectiveNeighbour(x,y,:),1,nAgeClasses))...
                    *transmissionConst* reshape(susceptible(x,y,:),nAgeClasses,1)/thisPop;
                rates(5) = rates(5) - transmissionRate(x,y) + new;
                transmissionRate(x,y)= new;
            else
                rates(5) = rates(5) - transmissionRate(x,y);
                transmissionRate(x,y) = 0;
            end

        case 6 % Birth
            [x,y] = ChoosePatch(birthRate);
            population(x,y,1) = population(x,y,1) + 1;
            susceptible(x,y,1) = susceptible(x,y,1) + 1;
            
            new = birthConst'*reshape(population(x,y,:),nAgeClasses,1) * max( [0, 1-sum(population(x,y,:))/K(x,y)]);
            rates(6) = rates(6) - birthRate(x,y) + new;
            birthRate(x,y) = new;
            new = sum((deathConst + migConst + 1) .* reshape(susceptible(x,y,:),nAgeClasses,1));
            rates(1) = rates(1) - DMRateS(x,y) + new;
            DMRateS(x,y) = new;
            
            thisPop = sum(population(x,y,:));
            if (thisPop > 0)
                new = (reshape(infective(x,y,:),1,nAgeClasses) + ...
                    globalTransmissionConst*reshape(infectiveNeighbour(x,y,:),1,nAgeClasses))...
                    *transmissionConst* reshape(susceptible(x,y,:),nAgeClasses,1)/thisPop;
                rates(5) = rates(5) - transmissionRate(x,y) + new;
                transmissionRate(x,y) = new;
                
            else
                rates(5) = rates(5) - transmissionRate(x,y);
                transmissionRate(x,y) = 0;
            end      
    end
    % Recording evolution of distance to disease front, prevalence(time
    % since arrival), quarterly field trips
    if ( (time-tMonitored) > 1/12)
        popNoAge = sum(population,3);
        infExp = sum((infective(:,:,:)+diagnosable(:,:,:)+exposed(:,:,:)),3);
        prevalence = 1 - sum(susceptible(:,:,3:end),3) ./ sum(population(:,:,3:end),3);
        
		nInf = [nInf; sum(sum(infExp))];
        timeRet = [timeRet; time];
        pop = [pop; sum(sum(popNoAge)) ];
        prev = [prev; mean(prevalence(~isnan(prevalence))) ];
        popInit = [popInit; popNoAge(xInit, yInit)];
        prevInit = [prevInit; prevalence(xInit, yInit)];
        occupied = [occupied; length(find(popNoAge > 0)) ];
        
        popSpatial(:,:,iRecord) = popNoAge;
        prevSpatial(:,:,iRecord) = prevalence;
        nInfSpatial(:,:,iRecord) = infExp;
        
        popAge(iRecord,:) = sum(sum(population));
        nInfAge (iRecord,:) = sum(sum((infective(:,:,:)+diagnosable(:,:,:)+exposed(:,:,:))));
        
        iRecord = iRecord + 1;

        % Stop if disease extinct or time past maximal time
        if ( sum(sum(infExp)) == 0 || sum(sum(popNoAge)) == 0)
            extinctionTime = time;
            break;
        end
        if (time > tMax)
            break;
        end
        
        tMonitored = time;
    end
end

prev(isnan(prev)) = 0;
prevInit(isnan(prevInit)) = 0;
prevSpatial(isnan(prevSpatial)) = 0;
