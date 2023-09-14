function flatStates = ecef2flat(refLL,states)

LLA = ecef2lla(states(7:9)','WGS84');
pos = lla2flat(LLA,refLL,0,0,'WGS84');
[u,v,w] = ecef2nedv(states(1),states(2),states(3),refLL(1),refLL(2),"degrees");

flatStates = [u;v;w;states(4);states(5);states(6);pos';states(10:end)];

end