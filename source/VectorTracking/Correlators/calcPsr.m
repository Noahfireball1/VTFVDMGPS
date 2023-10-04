function [psr,psrDot,unitVectors,activeSVs] = calcPsr(usrStates,svStates)
c = 299792458;
clkBias = usrStates(7);
clkDrift = usrStates(8);

dx = (svStates(1,:) - usrStates(1));
dy = (svStates(3,:) - usrStates(2));
dz = (svStates(5,:) - usrStates(3));

usrVel = usrStates(4:6);
svVel = [svStates(2,:);svStates(4,:);svStates(6,:)]';

range = sqrt(dx.^2 + dy.^2 + dz.^2);

psr = range;% + clkBias;

unitVectors = [dx./range; dy./range; dz./range];

for i = 1:length(psr)

    psrDot(i) =  unitVectors(:,i)'*(usrVel - svVel(i,:)'); %+ clkDrift;

end

% Discard any satellites with a negative elevation
LLA = ecef2lla(usrStates(1:3)','WGS84');
[~,el,~] = lookangles(LLA,[svStates(1,:);svStates(3,:);svStates(5,:)]');

activeSVs = el > 10;

end

