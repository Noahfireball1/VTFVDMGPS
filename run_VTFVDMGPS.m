%% Formatting
clc
clear
close all
format shortg

%% Adding Directories Based on User's Paths
projectRoot = fileparts(which(mfilename));
addpath(genpath(projectRoot))
dir.config = append(projectRoot,filesep,'config',filesep);
dir.source = append(projectRoot,filesep,'source',filesep);
dir.dataGPS = append(projectRoot,filesep,'data',filesep,'GPS',filesep);
dir.dataFVDM = append(projectRoot,filesep,'data',filesep,'FVDM',filesep);
dir.output = append(projectRoot,filesep,'output',filesep);

%% Initial Text to Console
utilities.printText.initial();

%% Selecting Configuration File
inputFile = uigetfile({'*.yaml'},'Select Input File',dir.config);
inputFilePath = append(dir.config,inputFile);

%% Configuring Vector Tracking Class
DCAM = VectorTracking(inputFilePath);
plotDCAM = VectorTrackingPlotting(inputFilePath);

%% Start Vector Tracking
utilities.printText(1)

for timeStep = 1:DCAM.timeStep

    DCAM.step(timeStep)

end