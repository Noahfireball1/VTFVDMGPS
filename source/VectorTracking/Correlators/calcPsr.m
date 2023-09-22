function [psr,psrDot,unitVectors,activeSVs] = calcPsr(usrStates,svStates)
c = 299792458;

if length(usrStates) < 14
    clkBias = 0;
    clkDrift = 0;
else
    clkBias = usrStates(13);
    clkDrift = usrStates(14);
end

dx = (svStates(1,:) - usrStates(7));
dy = (svStates(3,:) - usrStates(8));
dz = (svStates(5,:) - usrStates(9));

usrVel = usrStates(1:3);
svVel = [svStates(2,:);svStates(4,:);svStates(6,:)]';

range = sqrt(dx.^2 + dy.^2 + dz.^2);

psr = range + c*clkBias;

unitVectors = [dx./range; dy./range; dz./range];

for i = 1:length(psr)

    psrDot(i) =  unitVectors(:,i)'*(usrVel - svVel(i,:)') + c*clkDrift;

end

% Discard any satellites with a negative elevation
LLA = ecef2lla(usrStates(7:9)');
[~,el,~] = ecef2aer(svStates(1,:),svStates(3,:),svStates(5,:),LLA(1),LLA(2),LLA(3),wgs84Ellipsoid("meter"));

activeSVs = el > 20;

end

