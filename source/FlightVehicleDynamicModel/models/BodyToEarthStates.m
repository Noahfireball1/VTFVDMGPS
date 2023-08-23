function [x,y,z,u,v,w,clk,clkBias] = BodyToEarthStates(pos,vel,clk,clkBias,DCM_BL,refLL)

NED = DCM_BL*pos;
NEDv = DCM_BL*vel;

xNorth = NED(1);
yEast = NED(2);
zDown = NED(3);

uNorth = NEDv(1);
vEast = NEDv(2);
wDown = NEDv(3);

spheroid = wgs84Ellipsoid("meter");
[x,y,z] = ned2ecef(xNorth,yEast,zDown,refLL(1),refLL(2),0,spheroid);
[u,v,w] = ned2ecefv(uNorth,vEast,wDown,refLL(1),refLL(2));

end

