function [model,dir] = initializeSim(inputFilePath,dir)

% Initializes simulation variables, structures, etc...
printText(1)
config = ReadYaml(inputFilePath);
general = config.general;
aircraft = config.aircraft;
[~,configName,~] = fileparts(inputFilePath);

%% General
date = datetime(general.year,general.month,general.day);

if ~isempty(general.outputDir)
    dir.output = general.outputDir;
end

%% Model Related Properties

% Loading in required tables
printText(3)
load(append(dir.tables,"DA40.mat"));
printText(4)
load(append(dir.tables,"DA40ENGINE.mat"));
printText(5)
load(append(dir.tables,"DA40PF.mat"));
printText(6)
load(append(dir.tables,"DA40ST.mat"));
load(append(dir.tables,"DA40STGEOM.mat"));
load("DA40PROP.mat");
selWaypoints = load(sprintf('%s.mat',aircraft.waypoints));

% Pull Rinex for given date and time
rinex = GenerateEphemeris(date,dir);
rinex = string(rinex.rinexFilePath);
% Reading in Flight Vehicle's Initial State
u = aircraft.initialState.u;
v = aircraft.initialState.v;
w = aircraft.initialState.w;
p = aircraft.initialState.p;
q = aircraft.initialState.q;
r = aircraft.initialState.r;
x = aircraft.initialState.x;
y = aircraft.initialState.y;
z = aircraft.initialState.z;
phi = aircraft.initialState.phi;
theta = str2num(aircraft.initialState.theta);
psi = str2num(aircraft.initialState.psi);
clkBias = aircraft.initialState.clkBias;
clkDrift = aircraft.initialState.clkDrift;

printText(10)
fprintf('\t\t\t')
upd = textprogressbar(general.monteCarloRuns);
for i = 1:general.monteCarloRuns
    upd(i)
    rndSeedi = randi(1e6,1);
    model(i) = Simulink.SimulationInput('DA40_Flight_Model');
    model(i) = model(i).setVariable('rndSeed',rndSeedi);
    model(i) = model(i).setVariable('refLLA',selWaypoints.refLL);
    model(i) = model(i).setVariable('Start_Time',0);
    model(i) = model(i).setVariable('End_Time',general.duration);
    model(i) = model(i).setVariable('Time_Step',1/aircraft.frequency);
    model(i) = model(i).setVariable('Initial_X',[u;v;w;p;q;r;x;y;z;phi;theta;psi;clkBias;clkDrift]);
    model(i) = model(i).setVariable('rinexFilePath',rinex);
    model(i) = model(i).setVariable('dayofyear',day(date,"dayofyear"));
    model(i) = model(i).setVariable('waypoints',flipud(selWaypoints.waypoints));
    model(i) = model(i).setVariable('lookaheadDist',aircraft.lookaheadDistance);
    model(i) = model(i).setVariable('STGeometry',ST_Geometry);
    model(i) = model(i).setVariable('Vehicle',Vehicle);
    model(i) = model(i).setVariable('StNumbers',STNumbers);
    model(i) = model(i).setVariable('BSFC_LUT',BSFC_LUT);
    model(i) = model(i).setVariable('PowerFactor_LUT',PowerFactor_LUT);
    model(i) = model(i).setVariable('year',general.year);
    model(i) = model(i).setVariable('month',general.month);
    model(i) = model(i).setVariable('day',general.day);
    model(i) = model(i).setVariable('noisePowerForces',str2num(aircraft.noisePowerForces{1,1}));
    model(i) = model(i).setVariable('noiseTSForces',aircraft.noiseTSForces);
end

end

