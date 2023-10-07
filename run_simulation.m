%% Formatting
clc
clear
close all
format shortg
rng(1)
%% Adding Directories Based on User's Paths
projectRoot = fileparts(which(mfilename)); 
addpath(genpath(projectRoot))
dir.config = append(projectRoot,filesep,'config',filesep);
dir.src = append(projectRoot,filesep,'source',filesep);
dir.tables = append(dir.src,'tables',filesep);
dir.vt = append(dir.src,'VectorTracking',filesep);
dir.rinex = append(dir.vt,'Satellites',filesep,'rinex',filesep);
dir.svStates = append(dir.vt,'Satellites',filesep,'svStates',filesep);
dir.waypoints = append(projectRoot,filesep,'waypoints',filesep);
dir.output = append(projectRoot,filesep,'output',filesep);

%% Initial Text to Console
printText(8);

%% Selecting Configuration File
inputFile = uigetfile({'*.yaml'},'Select Input File',dir.config);
inputFilePath = append(dir.config,inputFile);

%% Initializing Simulation
[model,dir] = initializeSim(inputFilePath,dir);
%% Starting Simulation
printText(7);
run = sim(model);

save(append(dir.output,sprintf('%s_results.mat',inputFile(1:end-5))),"run")




