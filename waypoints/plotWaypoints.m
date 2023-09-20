%% Formatting
clc
clear
close all
format shortg

%% Load in Specified Set of Waypoints
load("SCurveFlightPath.mat")


LLA = flat2lla(waypoints,refLL,0,0,'WGS84');

figure('Position',[1000 200 900 800])
geoplot(LLA(:,1),LLA(:,2),'LineWidth',2,'Color',[0 0 0],LineStyle='--')
geobasemap satellite
ax = gca;
ax.FontSize = 16;