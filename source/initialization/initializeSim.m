% Initializes simulation variables, structures, etc...
config = ReadYaml(inputFilePath);
general = config.general;
aircraft = config.aircraft;

tableDir = append(dir.src,'tables',filesep);

date = datetime(general.year,general.month,general.day);
dayofyear = day(date,"dayofyear");

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


