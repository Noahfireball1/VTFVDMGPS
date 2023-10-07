function [model,dir] = initializeSim(inputFilePath,dir)

% Initializes simulation variables, structures, etc...
printText(1)
config = ReadYaml(inputFilePath);
general = config.general;
aircraft = config.aircraft;
signal = config.signal;
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
lat = aircraft.initialState.lat;
long = aircraft.initialState.long;
alt = aircraft.initialState.alt;
phi = str2num(aircraft.initialState.phi);
theta = str2num(aircraft.initialState.theta);
psi = str2num(aircraft.initialState.psi);


clkVar = formClkVariance(aircraft.clockType,1/aircraft.frequency);

clkNoise = sqrt(clkVar)*randn(2,1);
clkBias = clkNoise(1);
clkDrift = clkNoise(2);

initP = eye(14);
% initP(4:6,4:6) = 0;
% initP(10:12,10:12) = 0;
Q = blkdiag(diag([2500 2500 1000 1e-1 1e-1 1e-1 5e-10 5e-10 10000 0 0 0]),clkVar);
% initP = zeros(14);
% Q = zeros(14);

printText(10)
fprintf('\t\t\t')
upd = textprogressbar(general.monteCarloRuns);
for i = 1:general.monteCarloRuns
    upd(i)
    rndSeedi = randi(1e6,1);
    model(i) = Simulink.SimulationInput('DA40_Flight_Model');
    model(i) = model(i).setVariable('rndSeed',1);
    model(i) = model(i).setVariable('refLLA',selWaypoints.refLL);
    model(i) = model(i).setVariable('Start_Time',0);
    model(i) = model(i).setVariable('End_Time',general.duration);
    model(i) = model(i).setVariable('Time_Step',1/aircraft.frequency);
    model(i) = model(i).setVariable('Initial_X',[u;v;w;p;q;r;lat;long;alt;phi;theta;psi;clkBias;clkDrift]);
    model(i) = model(i).setVariable('rinexFilePath',rinex);
    model(i) = model(i).setVariable('dayofyear',day(date,"dayofyear"));
    model(i) = model(i).setVariable('waypoints',(selWaypoints.waypoints));
    model(i) = model(i).setVariable('lookaheadDist',aircraft.lookaheadDistance);
    model(i) = model(i).setVariable('STGeometry',ST_Geometry);
    model(i) = model(i).setVariable('Vehicle',Vehicle);
    model(i) = model(i).setVariable('MOI',Vehicle.MassProp.MOI);
    model(i) = model(i).setVariable('StNumbers',STNumbers);
    model(i) = model(i).setVariable('BSFC_LUT',BSFC_LUT);
    model(i) = model(i).setVariable('PowerFactor_LUT',PowerFactor_LUT);
    model(i) = model(i).setVariable('year',general.year);
    model(i) = model(i).setVariable('month',general.month);
    model(i) = model(i).setVariable('day',general.day);
    model(i) = model(i).setVariable('initCN0',10.^(signal.CN0/10).*ones(1,31));
    model(i) = model(i).setVariable('amplitude',signal.amplitude);
    model(i) = model(i).setVariable('variance',str2num(cell2mat(aircraft.noiseVariance)));
    model(i) = model(i).setVariable('clkVar',clkVar);
    model(i) = model(i).setVariable('initP',initP);
    model(i) = model(i).setVariable('Q',Q);
end

end

