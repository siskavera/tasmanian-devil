function [newX,newY] = ChoosePatchForMigrant(isSea, x, y)
% CHOOSEPATCHFORMIGRANT Choose patch for a migrant
%   ChoosePatchForMigrant(isSea, x, y) Returns indices for the new patch of
%   a migrant devil. New patch chosen uniformly at random, from valid (not isSea)
%   patches out of the neighbouring 4 patches
newX = x;
newY = y;

while (true)
    switch ceil(4*rand)
        case 1
            newX = x - 1;
        case 2
            newX = x + 1;
        case 3
            newY = y - 1;
        case 4
            newY = y + 1;
    end
    if (~isSea(newX,newY))
        break;
    end
    newX = x;
    newY = y;
end