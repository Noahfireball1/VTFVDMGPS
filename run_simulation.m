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
printText(8);

%% Selecting Configuration File
inputFile = uigetfile({'*.yaml'},'Select Input File',dir.config);
inputFilePath = append(dir.config,inputFile);

%% Initializing Simulation
model = initializeSim(inputFilePath,dir);
%% Starting Simulation
printText(7);

run = parsim(in);





