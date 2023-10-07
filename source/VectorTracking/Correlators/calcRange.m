function [range,rangeRate] = calcRange(usrStates,svStates)

usrVel = [usrStates(2);usrStates(4);usrStates(6)];
svVel = [svStates(2,:);svStates(4,:);svStates(6,:)];

range = sqrt((svStates(1,:) - usrStates(1)).^2 ...
    + (svStates(3,:) - usrStates(3)).^2 ...
    + (svStates(1,:) - usrStates(5)).^2);

unitVectors = [(svStates(1,:) - usrStates(1))./range;...
    (svStates(3,:) - usrStates(3))./range;...
    (svStates(5,:) - usrStates(5))./range];

for i = 1:length(range)

    ux = unitVectors(1,i);
    uy = unitVectors(2,i);
    uz = unitVectors(3,i);

    rangeRate(i) =  -ux*(svVel(1,i) - usrVel(1)) + ...
        -uy*(svVel(2,i) - usrVel(2)) + ...
        -uz*(svVel(3,i) - usrVel(3));

end