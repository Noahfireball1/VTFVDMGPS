function [psr,psrDot,unitVectors,activeSVs] = calcPsr(usrStates,svStates)

clkBias = usrStates(7);
clkDrift = usrStates(8);

usrVel = [usrStates(2);usrStates(4);usrStates(6)];
svVel = [svStates(2,:);svStates(4,:);svStates(6,:)];

range = sqrt((svStates(1,:) - usrStates(1)).^2 ...
    + (svStates(3,:) - usrStates(3)).^2 ...
    + (svStates(1,:) - usrStates(5)).^2) - clkBias;

psr = range + clkBias;

unitVectors = [(svStates(1,:) - usrStates(1))./range;...
    (svStates(3,:) - usrStates(3))./range;...
    (svStates(5,:) - usrStates(5))./range];

for i = 1:length(psr)

    ux = unitVectors(1,i);
    uy = unitVectors(2,i);
    uz = unitVectors(3,i);

    psrDot(i) =  ux*(svVel(1,i) - usrVel(1)) + ...
        uy*(svVel(2,i) - usrVel(2)) + ...
        uz*(svVel(3,i) - usrVel(3)) + clkDrift;

end

% Discard any satellites with a negative elevation
LLA = ecef2lla([usrStates(1) usrStates(3) usrStates(5)],'WGS84');
[~,el,~] = lookangles(LLA,[svStates(1,:);svStates(3,:);svStates(5,:)]');

activeSVs = el > 10;

end