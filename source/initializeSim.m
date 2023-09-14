function model = initializeSim(inputFilePath,dir)

% Initializes simulation variables, structures, etc...
config = ReadYaml(inputFilePath);
general = config.general;
aircraft = config.aircraft;

tableDir = append(dir.src,'tables',filesep);

date = datetime(general.year,general.month,general.day);
dayofyear = day(date,"dayofyear");


% Generating Satellite States
satStates = genSatellitesStates(End_Time,date,dir);

for i = 1:1
    rndSeedi = randi(1e6,1);
    in(i) = Simulink.SimulationInput('DA40_Flight_Model');
    in(i) = in(i).setVariable('rndSeed',rndSeedi);
    in(i) = in(i).setVariable('refLLA',refLLA);
    in(i) = in(i).setVariable('Start_Time',Start_Time);
    in(i) = in(i).setVariable('End_Time',End_Time);
    in(i) = in(i).setVariable('Time_Step',Time_Step);
    in(i) = in(i).setVariable('Initial_X',Initial_X);
    in(i) = in(i).setVariable('satStates',satStates);
    in(i) = in(i).setVariable('dayofyear',dayofyear);
    in(i) = in(i).setVariable('waypoints',waypoints);
    in(i) = in(i).setVariable('lookaheadDist',lookaheadDist);
    in(i) = in(i).setVariable('STGeometry',STGeometry);
    in(i) = in(i).setVariable('Vehicle',Vehicle);
    in(i) = in(i).setVariable('StNumbers',STNumbers);
    in(i) = in(i).setVariable('BSFC_LUT',BSFC_LUT);
    in(i) = in(i).setVariable('PowerFactor_LUT',PowerFactor_LUT);
end

Start_Time = 0;
End_Time = aircraft.time;
Time_Step = 1/aircraft.frequency;

if ~isempty(general.outputDir)
    outputDir = general.outputDir;
end

if config.general.verbose
    printText(1);
else
    printText(2);
end

if config.general.verbose
    printText(3);
end
Vehicle = load(append(tableDir,"DA40.mat"));
BSFC_LUT = load(append(tableDir,"DA40ENGINE.mat"));
PowerFactor_LUT = load(append(tableDir,"DA40PF.mat"));

Vehicle = Vehicle.Vehicle;
BSFC_LUT = BSFC_LUT.BSFC_LUT;
PowerFactor_LUT = PowerFactor_LUT.PowerFactor_LUT;

if config.general.verbose
    printText(4);
end

STNumbers = load(append(tableDir,"DA40ST.mat"));
STGeometry = load(append(tableDir,"DA40STGEOM.mat"));

STNumbers = STNumbers.STNumbers;
STGeometry = STGeometry.ST_Geometry;

if config.general.verbose
    printText(5);
end

Prop = load("DA40PROP.mat");
Prop = Prop.Prop;

selWaypoints = load(sprintf('%s.mat',config.aircraft.waypoints));
waypoints = selWaypoints.waypoints;
lookaheadDist = aircraft.lookaheadDistance;


refLLA = selWaypoints.refLL;
Initial_X = [aircraft.initialState.u;...
    aircraft.initialState.v;...
    aircraft.initialState.w;...
    aircraft.initialState.p;...
    aircraft.initialState.q;...
    aircraft.initialState.r;...
    aircraft.initialState.x;...
    aircraft.initialState.y;...
    aircraft.initialState.z;...
    aircraft.initialState.phi;...
    str2num(aircraft.initialState.theta);...
    str2num(aircraft.initialState.psi);...
    0;...
    0];

if config.general.verbose
    printText(6);
end
end

