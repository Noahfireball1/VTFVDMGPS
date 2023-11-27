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
dir.tables = append(dir.src,'tables',filesep);
dir.vt = append(dir.src,'VectorTracking',filesep);
dir.rinex = append(dir.vt,'Satellites',filesep,'rinex',filesep);
dir.svStates = append(dir.vt,'Satellites',filesep,'svStates',filesep);
dir.waypoints = append(projectRoot,filesep,'waypoints',filesep);
dir.output = append(projectRoot,filesep,'output',filesep);

%% Initial Text to Console
printText(8);

%% 45 db MC
% Initializing Simulation
[model,dir] = initializeSim('C:\Users\noahm\graduate\thesis\VTFVDMGPS\config\MC_100_45db.yaml',dir);
% Starting Simulation
printText(7);
run = parsim(model);

save(append(dir.output,'MC_100_45db_results.mat'),"run","-v7.3")

clear run

%% 35 db MC
% Initializing Simulation
[model,dir] = initializeSim('C:\Users\noahm\graduate\thesis\VTFVDMGPS\config\MC_100_35db.yaml',dir);
% Starting Simulation
printText(7);
run = parsim(model);

save(append(dir.output,'MC_100_35db_results.mat'),"run","-v7.3")

clear run
%% 25 db MC
% Initializing Simulation
[model,dir] = initializeSim('C:\Users\noahm\graduate\thesis\VTFVDMGPS\config\MC_100_25db.yaml',dir);
% Starting Simulation
printText(7);
run = parsim(model);

save(append(dir.output,'MC_100_25db_results.mat'),"run","-v7.3")

clear run
%% 20 db MC
% Initializing Simulation
[model,dir] = initializeSim('C:\Users\noahm\graduate\thesis\VTFVDMGPS\config\MC_100_20db.yaml',dir);
% Starting Simulation
printText(7);
run = parsim(model);

save(append(dir.output,'MC_100_20db_results.mat'),"run","-v7.3")

clear run
%% 15 db MC
% Initializing Simulation
[model,dir] = initializeSim('C:\Users\noahm\graduate\thesis\VTFVDMGPS\config\MC_100_15db.yaml',dir);
% Starting Simulation
printText(7);
run = parsim(model);

save(append(dir.output,'MC_100_15db_results.mat'),"run","-v7.3")

clear run
%% 10 db MC
% Initializing Simulation
[model,dir] = initializeSim('C:\Users\noahm\graduate\thesis\VTFVDMGPS\config\MC_100_10db.yaml',dir);
% Starting Simulation
printText(7);
run = parsim(model);

save(append(dir.output,'MC_100_10db_results.mat'),"run","-v7.3")

clear run
%% 5 db MC
% Initializing Simulation
[model,dir] = initializeSim('C:\Users\noahm\graduate\thesis\VTFVDMGPS\config\MC_100_5db.yaml',dir);
% Starting Simulation
printText(7);
run = parsim(model);

save(append(dir.output,'MC_100_5db_results.mat'),"run","-v7.3")

clear run
%% 2 db MC
% Initializing Simulation
[model,dir] = initializeSim('C:\Users\noahm\graduate\thesis\VTFVDMGPS\config\MC_100_2db.yaml',dir);
% Starting Simulation
printText(7);
run = parsim(model);

save(append(dir.output,'MC_100_2db_results.mat'),"run","-v7.3")

clear run