function [range,rangeRate] = calcRange(usrStates,svStates)

dx = (svStates(1,:) - usrStates(1));
dy = (svStates(3,:) - usrStates(2));
dz = (svStates(5,:) - usrStates(3));

usrVel = usrStates(4:6)';
svVel = [svStates(2,:);svStates(4,:);svStates(6,:)];

range = sqrt(dx.^2 + dy.^2 + dz.^2);

unitVectors = [dx./range; dy./range; dz./range];

for i = 1:length(range)

    ux = unitVectors(1,i);
    uy = unitVectors(2,i);
    uz = unitVectors(3,i);

    rangeRate(i) =  ux*(usrVel(1) - svVel(1,i)) + ...
        uy*(usrVel(2) - svVel(2,i)) + ...
        uz*(usrVel(3) - svVel(3,i));

end