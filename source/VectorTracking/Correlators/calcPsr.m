function [psr,carrFreq,unitVectors] = calcPsr(usrStates,svStates)

dx = (svStates(1,:) - usrStates(7));
dy = (svStates(3,:) - usrStates(8));
dz = (svStates(5,:) - usrStates(9));
usrVel = usrStates(1:3);
svVel = [svStates(2,:);svStates(4,:);svStates(6,:)]';

range = sqrt(dx.^2 + dy.^2 + dz.^2);

psr = range;
unitVectors = [dx./range; dy./range; dz./range];

for i = 1:length(psr)

    carrFreq(i) =  1575.42e6 - (299792458/1575.42e6)*(svVel(i,:)*unitVectors(:,i)) + (299792458/1575.42e6)*(usrVel'*unitVectors(:,i));

end

% Discard any satellites with a negative elevation
LLA = ecef2lla(usrStates(7:9)');
[~,el,~] = ecef2aer(svStates(1,:),svStates(3,:),svStates(5,:),LLA(1),LLA(2),LLA(3),wgs84Ellipsoid("meter"));

psr = psr(1,el > 0);
carrFreq = carrFreq(1,el > 0);
unitVectors = unitVectors(:,el > 0);

end

