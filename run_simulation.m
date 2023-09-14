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

run = sim("DA40_Flight_Model.slx");

estLLA = flat2lla(run.estimatedStates(:,7:9),refLLA,0,0,'WGS84');
refLLA = flat2lla(run.rcvrStates(:,7:9),refLLA,0,0,'WGS84');


%% Plotting
figure
geoplot(refLLA(:,1),refLLA(:,2))
hold on 
geoplot(estLLA(:,1),estLLA(:,2))





