% Initializes simulation variables, structures, etc...
config = ReadYaml(inputFilePath);
general = config.general;
aircraft = config.aircraft;
noise = config.noise;

tableDir = append(srcDir,'tables',filesep);

date = datetime(general.year,general.month,general.day);
dayofyear = day(date,"dayofyear");

Start_Time = 0;
End_Time = aircraft.time;
Time_Step = 1/aircraft.frequency;
imuTimeStep = 1/noise.IMUFrequency;

numRuns = general.numMonteCarloRuns;
saveResults = general.save;

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
    str2num(aircraft.initialState.psi)];

if config.general.verbose
    printText(6);
end

wind = config.noise.FVDM.wind;
gravity = config.noise.FVDM.gravity;
gravity.force = str2num(cell2mat(gravity.force));
gravity.moment = str2num(cell2mat(gravity.moment));

aero = config.noise.FVDM.aero;
aero.force = str2num(cell2mat(aero.force));
aero.moment = str2num(cell2mat(aero.moment));

engine = config.noise.FVDM.engine;
engine.force = str2num(cell2mat(engine.force));
engine.moment = str2num(cell2mat(engine.moment));

IMU1 = config.noise.IMU1;
IMU1.accelerometer.accelerationRandomWalk = str2num(cell2mat(IMU1.accelerometer.accelerationRandomWalk));
IMU1.accelerometer.axesSkew = str2num(cell2mat(IMU1.accelerometer.axesSkew));
IMU1.accelerometer.biasInstability = str2num(cell2mat(IMU1.accelerometer.biasInstability));
IMU1.accelerometer.constantOffsetBias = str2num(cell2mat(IMU1.accelerometer.constantOffsetBias));
IMU1.accelerometer.temperatureBias = str2num(cell2mat(IMU1.accelerometer.temperatureBias));
IMU1.accelerometer.temperatureScaleFactor = str2num(cell2mat(IMU1.accelerometer.temperatureScaleFactor));
IMU1.accelerometer.velocityRandomWalk = str2num(cell2mat(IMU1.accelerometer.velocityRandomWalk));
IMU1.gyroscope.accelerationBias = str2num(cell2mat(IMU1.gyroscope.accelerationBias));
IMU1.gyroscope.angleRandomWalk = str2num(cell2mat(IMU1.gyroscope.angleRandomWalk));
IMU1.gyroscope.axesSkew = str2num(cell2mat(IMU1.gyroscope.axesSkew));
IMU1.gyroscope.biasInstability = str2num(cell2mat(IMU1.gyroscope.biasInstability));
IMU1.gyroscope.constantOffsetBias = str2num(cell2mat(IMU1.gyroscope.constantOffsetBias));
IMU1.gyroscope.gyroscopeRandomWalk = str2num(cell2mat(IMU1.gyroscope.gyroscopeRandomWalk));
IMU1.gyroscope.temperatureBias = str2num(cell2mat(IMU1.gyroscope.temperatureBias));
IMU1.gyroscope.temperatureScaleFactor = str2num(cell2mat(IMU1.gyroscope.temperatureScaleFactor));

IMU2 = config.noise.IMU2;
IMU2.accelerometer.accelerationRandomWalk = str2num(cell2mat(IMU2.accelerometer.accelerationRandomWalk));
IMU2.accelerometer.axesSkew = str2num(cell2mat(IMU2.accelerometer.axesSkew));
IMU2.accelerometer.biasInstability = str2num(cell2mat(IMU2.accelerometer.biasInstability));
IMU2.accelerometer.constantOffsetBias = str2num(cell2mat(IMU2.accelerometer.constantOffsetBias));
IMU2.accelerometer.temperatureBias = str2num(cell2mat(IMU2.accelerometer.temperatureBias));
IMU2.accelerometer.temperatureScaleFactor = str2num(cell2mat(IMU2.accelerometer.temperatureScaleFactor));
IMU2.accelerometer.velocityRandomWalk = str2num(cell2mat(IMU2.accelerometer.velocityRandomWalk));
IMU2.gyroscope.accelerationBias = str2num(cell2mat(IMU2.gyroscope.accelerationBias));
IMU2.gyroscope.angleRandomWalk = str2num(cell2mat(IMU2.gyroscope.angleRandomWalk));
IMU2.gyroscope.axesSkew = str2num(cell2mat(IMU2.gyroscope.axesSkew));
IMU2.gyroscope.biasInstability = str2num(cell2mat(IMU2.gyroscope.biasInstability));
IMU2.gyroscope.constantOffsetBias = str2num(cell2mat(IMU2.gyroscope.constantOffsetBias));
IMU2.gyroscope.gyroscopeRandomWalk = str2num(cell2mat(IMU2.gyroscope.gyroscopeRandomWalk));
IMU2.gyroscope.temperatureBias = str2num(cell2mat(IMU2.gyroscope.temperatureBias));
IMU2.gyroscope.temperatureScaleFactor = str2num(cell2mat(IMU2.gyroscope.temperatureScaleFactor));

IMU3 = config.noise.IMU3;
IMU3.accelerometer.accelerationRandomWalk = str2num(cell2mat(IMU3.accelerometer.accelerationRandomWalk));
IMU3.accelerometer.axesSkew = str2num(cell2mat(IMU3.accelerometer.axesSkew));
IMU3.accelerometer.biasInstability = str2num(cell2mat(IMU3.accelerometer.biasInstability));
IMU3.accelerometer.constantOffsetBias = str2num(cell2mat(IMU3.accelerometer.constantOffsetBias));
IMU3.accelerometer.temperatureBias = str2num(cell2mat(IMU3.accelerometer.temperatureBias));
IMU3.accelerometer.temperatureScaleFactor = str2num(cell2mat(IMU3.accelerometer.temperatureScaleFactor));
IMU3.accelerometer.velocityRandomWalk = str2num(cell2mat(IMU3.accelerometer.velocityRandomWalk));
IMU3.gyroscope.accelerationBias = str2num(cell2mat(IMU3.gyroscope.accelerationBias));
IMU3.gyroscope.angleRandomWalk = str2num(cell2mat(IMU3.gyroscope.angleRandomWalk));
IMU3.gyroscope.axesSkew = str2num(cell2mat(IMU3.gyroscope.axesSkew));
IMU3.gyroscope.biasInstability = str2num(cell2mat(IMU3.gyroscope.biasInstability));
IMU3.gyroscope.constantOffsetBias = str2num(cell2mat(IMU3.gyroscope.constantOffsetBias));
IMU3.gyroscope.gyroscopeRandomWalk = str2num(cell2mat(IMU3.gyroscope.gyroscopeRandomWalk));
IMU3.gyroscope.temperatureBias = str2num(cell2mat(IMU3.gyroscope.temperatureBias));
IMU3.gyroscope.temperatureScaleFactor = str2num(cell2mat(IMU3.gyroscope.temperatureScaleFactor));


function printText(option)

programStyle = 'green';
textStyle = '*green';

switch option
    case 1
        cprintf(programStyle,'[DCAM]\t')
        cprintf(textStyle, 'Loading selected configuration file\n')
        cprintf(textStyle, '\t\tRunning in Verbose mode...\n\n')

    case 2
        cprintf(programStyle,'[DCAM]\t')
        cprintf(textStyle, 'Loading selected configuration file\n')
        cprintf(textStyle, '\t\tRunning in Silent mode...\n\n')
    case (3)
        cprintf(programStyle,'[DCAM]\t')
        cprintf(textStyle, 'Loading Diamond DA-40 aircraft properties...\n')
    case(4)
        cprintf(programStyle,'[DCAM]\t')
        cprintf(textStyle, 'Loading Diamond DA-40 aerodynamic properties...\n')
    case(5)
        cprintf(programStyle,'[DCAM]\t')
        cprintf(textStyle, 'Loading 3-blade propeller properties...\n')
    case(6)
        cprintf(programStyle,'[DCAM]\t')
        cprintf(textStyle, 'Loading noise models for FVDM and selected IMUs...\n')

    otherwise
end

end


