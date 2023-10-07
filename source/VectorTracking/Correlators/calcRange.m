function [range,rangeRate,unitVectors,activeSVs] = calcRange(usrPos,usrVel,svPos,svVel)
rangeRate = zeros(31,1);

range = sqrt((svPos(:,1) - usrPos(1)).^2 ...
    + (svPos(:,2) - usrPos(2)).^2 ...
    + (svPos(:,3) - usrPos(3)).^2);


unitVectors = [(svPos(:,1) - usrPos(1))./range,...
    (svPos(:,2) - usrPos(2))./range,...
    (svPos(:,3) - usrPos(3))./range];

for i = 1:31

    ux = unitVectors(i,1);
    uy = unitVectors(i,2);
    uz = unitVectors(i,3);

    rangeRate(i) =  ux*(usrVel(1) - svVel(i,1)) + ...
        uy*(usrVel(2) - svVel(i,2)) + ...
        uz*(usrVel(3) - svVel(i,3));

end

% Discard any satellites with a negative elevation
LLA = ecef2lla(usrPos,'WGS84');
[~,el,~] = lookangles(LLA,svPos);

activeSVs = el > 10;