%% Formatting
clc
clear
close all
format shortg
rng(1,"v5normal")

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

%% Load in Simulation Results
simulation = load("example_simulation.mat");

f_ib_b = simulation.run.Forces;
m_ib_b = simulation.run.Moments;

rcvrStates = simulation.run.truthStates';

%% Setting Configuration
year = 2022;
month = 12;
day = 22;
date = datetime(year,month,day);
duration = 30;
frequency = 200;
clockType = 'OCXO';
initialState = [32.609856 0 -85.480782 0 200 0 0 0];
CN0 = 45;

% Pull Rinex for given date and time
rinex = GenerateEphemeris(date,dir);
rinexFilePath = string(rinex.rinexFilePath);

%% Begin Vector Tracking
xMinus(:,1) = initialState;
clkVar = formClkVariance(clockType,1/frequency);
Q = blkdiag(diag([1.5 0 1.5 0 3.0 0]),clkVar);
pMinus(:,:,1) = Q;
timeArray = 0:1/frequency:duration;

signal.initCN0 = ones(1,31).*10^(CN0/10);
signal.phase = unifrnd(0,1,1,31);
signal.noise = (1./(2*(1/50)*signal.initCN0));
signal.amplitude = 4*signal.noise;
signal.CN0 = signal.initCN0;

for timeIdx = 2:length(timeArray)

    % Time Update
    [F,G] = formFG(xMinus(:,timeIdx-1),1/frequency);
    pMinus(:,:,timeIdx) = F*pMinus(:,:,timeIdx-1)*F' + F*Q*F'*(1/frequency);
    xMinus(:,timeIdx) = predictStates(xMinus(:,timeIdx-1),f_ib_b(timeIdx,:)',m_ib_b(timeIdx,:)',1/frequency,clkVar);

    % Measurement Update
    if mod(timeArray(timeIdx),1/50) == 0

        % Set Current Receiver States
        rcvr = rcvrStates(:,timeIdx);

        % Set Current Predicted States
        xTmp = xMinus(:,timeIdx);

        % Convert Receiver States to ECEF
        refStates(1:3,1) = lla2ecef(rcvr(7:9)'.*(180/pi),'WGS84')';
        refStates(4:6,1) = ned2ecefv(rcvr(1),rcvr(2),rcvr(3),rcvr(7),rcvr(8),'radians');
        refStates(7,1) = 0;
        refStates(8,1) = 0;

        % Convert Predicted States to ECEF
        estStates(1:3,1) = lla2ecef(xTmp(7:9)'.*(180/pi),'WGS84');
        estStates(4:6,1) = ned2ecefv(xTmp(1),xTmp(2),xTmp(3),xTmp(7),xTmp(8),'radians');
        estStates(7,1) = xTmp(13);
        estStates(8,1) = xTmp(14);

        % Pull Current Satellite States
        svStates = genSatellitesStates(timeArray(timeIdx),year,month,day,rinexFilePath);

        % Calculate Reciever PSR, Carrier Frequency, Unit Vectors
        [psr,carr,~,activeSVs] = calcPsr(refStates,svStates);
        refPsr = psr(activeSVs);
        refCarr = carr(activeSVs);

        % Calculate Predicted PSR, Carrier Frequency
        [psr,carr,unitVectors,~] = calcPsr(estStates,svStates);
        estPsr = psr(activeSVs);
        estCarr = carr(activeSVs);
        unitVectors = unitVectors(:,activeSVs);

        % Generate Correlator Residuals
        [resPsr,resCarr,varPsr,varCarr,signal] = genCorrelatorResiduals(refPsr,estPsr,refCarr,estCarr,signal,activeSVs);

        % Form Z Array
        Z = formZ(resPsr,resCarr);

        % Form H Matrix
        H = formH(unitVectors,xTmp);

        % Form R Matrix
        R = formR(varPsr,varCarr);

        % Update Kalman Gain
        L = calcL(H,pMinus(:,:,timeIdx),R);

        % Update Estimated States
        xPlus = xTmp - L*Z;

        % Update Estimated Covariance
        pPlus = updateCovariance(pMinus(:,:,timeIdx),L,H,R);

        xMinus(:,timeIdx) = xPlus;
        pMinus(:,:,timeIdx) = pPlus;

    end

end
figure
geoplot(xMinus(7,:).*(180/pi),xMinus(8,:).*(180/pi));
hold on
geoplot(rcvrStates(7,:).*(180/pi),rcvrStates(8,:).*(180/pi))