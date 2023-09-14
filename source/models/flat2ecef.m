function ecefStates = flat2ecef(refLL,states)

LLA = flat2lla(states(7:9)',refLL,0,0,'WGS84');
pos = lla2ecef(LLA,'WGS84');
[u,v,w] = ned2ecefv(states(1),states(2),states(3),refLL(1),refLL(2),"degrees");

ecefStates = [u;v;w;states(4);states(5);states(6);pos';states(10:end)];

end