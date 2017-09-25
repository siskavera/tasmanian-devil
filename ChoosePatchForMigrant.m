function [newX,newY] = ChoosePatchForMigrant(isSea, x, y)

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