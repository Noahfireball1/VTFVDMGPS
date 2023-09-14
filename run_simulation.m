%% Formatting
clc
clear
close all
format shortg

%% Adding Directories Based on User's Paths
projectRoot = fileparts(which(mfilename));
addpath(genpath(projectRoot))
dir.config = append(projectRoot,filesep,'config',filesep);
dir.src = append(projectRoot,filesep,'source',filesep);
dir.rinex = append(dir.src,'Satellites',filesep,'rinex',filesep);
dir.waypoints = append(projectRoot,filesep,'waypoints',filesep);
dir.output = append(projectRoot,filesep,'output',filesep);

%% Initial Text to Console
initialText;

%% Selecting Configuration File
inputFile = uigetfile({'*.yaml'},'Select Input File',dir.config);
inputFilePath = append(dir.config,inputFile);

%% Initializing Simulation
initializeSim;
%% Starting Simulation
simText;

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

run = parsim(in);





