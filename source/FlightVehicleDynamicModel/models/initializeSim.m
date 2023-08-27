% Initializes simulation variables, structures, etc...

config = ReadYaml(inputFilePath);
date = datetime(config.general.year,config.general.month,config.general.day);

Vehicle = load("DA40.mat");
Vehicle = Vehicle.Vehicle;
BSFC_LUT = load("DA40ENGINE.mat");
BSFC_LUT = BSFC_LUT.BSFC_LUT;
PowerFactor_LUT = load("DA40PF.mat");
PowerFactor_LUT = PowerFactor_LUT.PowerFactor_LUT;

STGeometry = load("DA40STGEOM.mat");
STGeometry = STGeometry.ST_Geometry;

Prop = load("DA40PROP.mat");
Prop = Prop.Prop;


startTime = 0;
endTime = config.aircraft.time;
timeStep = 1/config.aircraft.frequency;
selWaypoints = load(sprintf('%s.mat',config.aircraft.waypoints));
waypoints = selWaypoints.waypoints;
dayofyear = day(date,'dayofyear');
lookaheadDist = config.aircraft.lookaheadDistance;
refLLA = selWaypoints.refLL;


Initial_X = [config.aircraft.initialState.u;...
    config.aircraft.initialState.v;...
    config.aircraft.initialState.w;...
    config.aircraft.initialState.p;...
    config.aircraft.initialState.q;...
    config.aircraft.initialState.r;...
    config.aircraft.initialState.x;...
    config.aircraft.initialState.y;...
    config.aircraft.initialState.z;...
    config.aircraft.initialState.phi;...
    str2num(config.aircraft.initialState.theta);...
    str2num(config.aircraft.initialState.psi)];



