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
clkBias = aircraft.initialState.clkBias;
clkDrift = aircraft.initialState.clkDrift;

clkVar = formClkVariance(aircraft.clockType,1/aircraft.frequency);

B = zeros(12,6);
B(1:6,1:6) = eye(6);

initP = [      0.55192     0.026448    -0.020345    0.0025595      0.13879       0.1032   1.0106e-09   -4.271e-10   0.00074523     0.017942    -0.020902     0.043223  -7.3945e-05  -0.00017436
     0.026448      0.54021    -0.022944     0.035876      0.10404     0.082959   1.5543e-10  -2.0828e-09    0.0016732      0.01916     -0.01662     0.026338   -0.0001609  -0.00042476
    0.020345    0.022944     0.029221    -0.040165     -0.14457     -0.11247  -1.4504e-10   6.1264e-10   -0.0028468    -0.025328     0.022434     -0.03881   0.00044076     0.001217
    0.0025595     0.035876    -0.040165      0.42856     0.043181     0.065733  -2.1312e-11  -6.0985e-11   4.8295e-05     0.063898   -0.0036588    0.0084604  -1.6746e-07    1.051e-06
      0.13879      0.10404     -0.14457     0.043181       1.3328      0.83955   1.8622e-10  -8.0159e-10  -0.00043114      0.16035     -0.15682      0.32432  -7.3251e-05   6.9882e-05
       0.1032     0.082959     -0.11247     0.065733      0.83955      0.89899   1.5447e-10   -6.078e-10  -0.00024656      0.11067     -0.15405      0.21521  -4.6269e-05   3.7917e-05
   1.0106e-09   1.5543e-10  -1.4504e-10  -2.1312e-11   1.8622e-10   1.5447e-10   2.2324e-14  -3.9746e-14   1.2593e-07   1.7059e-11  -3.5827e-11   1.3213e-10  -2.7252e-08  -2.5964e-11
   -4.271e-10  -2.0828e-09   6.1264e-10  -6.0985e-11  -8.0159e-10   -6.078e-10  -3.9746e-14   1.4934e-13  -5.0775e-07  -1.0916e-10   1.3369e-10  -3.6118e-10   1.1985e-07   1.1317e-10
   0.00074523    0.0016732   -0.0028468   4.8295e-05  -0.00043114  -0.00024656   1.2593e-07  -5.0775e-07       2.9958  -3.8875e-05   1.9154e-07  -0.00011459      -1.1422   -0.0010823
     0.017942      0.01916    -0.025328     0.063898      0.16035      0.11067   1.7059e-11  -1.0916e-10  -3.8874e-05       0.0325    -0.021171     0.043371  -7.8143e-06   8.2961e-06
    -0.020902     -0.01662     0.022434   -0.0036588     -0.15682     -0.15405  -3.5827e-11   1.3369e-10   1.9114e-07    -0.021171     0.030993    -0.042872   6.6412e-06   2.5204e-07
     0.043223     0.026338     -0.03881    0.0084604      0.32432      0.21521   1.3213e-10  -3.6118e-10  -0.00011459     0.043371    -0.042872      0.10512  -1.9131e-05   1.8813e-05
  -7.3945e-05   -0.0001609   0.00044076  -1.6736e-07  -7.3248e-05  -4.6267e-05  -2.7252e-08   1.1985e-07      -1.1422   -7.814e-06   6.6409e-06  -1.9131e-05      2e6   2e3
  -0.00017436  -0.00042476     0.001217    1.051e-06   6.9882e-05   3.7917e-05  -2.5964e-11   1.1317e-10   -0.0010823   8.2961e-06   2.5204e-07   1.8813e-05   2e3   2e6];


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
end

end

