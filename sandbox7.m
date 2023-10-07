%% Formatting
clc
clear
close all
format shortg
rng(1)
rinexFilePath = "C:\Users\nsm0014\graduate\thesis\VTFVDMGPS\source\VectorTracking\Satellites\rinex\AMC400USA_R_20223560000_01D_GN.rnx";
%% Static Vector Tracking Analysis
AuburnLLA = [32.6099,-85.4808,200];
AuburnECEF = lla2ecef(AuburnLLA,'WGS84');
x_m(:,1) = [AuburnECEF(1);0;AuburnECEF(2);0;AuburnECEF(3);0;randn(1)*1e-10;randn(1)*1e-12];
timeStep = 1/50;
time = 0:timeStep:50-timeStep;
k = [0 1;0 0];
A = blkdiag(k,k,k,k);
Phi = expm(A*timeStep);

c = 299792458;
h_0 = 2e-25;
h_2 = 6e-25;
sigmaX = 1;
sigmaY = 1;
sigmaZ = 1;
S_f = c^2*h_0/2;
S_g = c^2*2*pi^2*h_2;
Q1 = [timeStep^3/3 timeStep^2/2; timeStep^2/2 timeStep]*sigmaX;
Q2 = [timeStep^3/3 timeStep^2/2; timeStep^2/2 timeStep]*sigmaY;
Q3 = [timeStep^3/3 timeStep^2/2; timeStep^2/2 timeStep]*sigmaZ;
Q4 = [S_f*timeStep + S_g*timeStep^3/3 S_g*timeStep^2/2;S_g*timeStep^2/2 S_g*timeStep];
Qd = blkdiag(Q1,Q2,Q3,Q4);

p_m(:,:,1) = eye(8);

truthStates(:,1) = x_m(:,1);

initCN0 = 10^(50/10).*ones(31,1);
oldAmplitude = ones(31,1);
oldNoise = 1./(2*0.02*(initCN0));
oldPhase = unifrnd(0,1,31,1);
oldCN0 = initCN0;



for timeIdx = 2:length(time)

    satStates = genSatellitesStates(time(timeIdx),2022,12,22,rinexFilePath);
    svPos = [satStates(1,:)' satStates(3,:)' satStates(5,:)'];
    svVel = [satStates(2,:)' satStates(4,:)' satStates(6,:)'];
    %% Truth State Propagation
    truthStates(:,timeIdx) = Phi*truthStates(:,timeIdx-1) + sqrt(Qd)*randn(8,1);

    %% Time Update
    x_m(:,timeIdx) = Phi*x_m(:,timeIdx-1);
    p_m(:,:,timeIdx) = Phi*p_m(:,:,timeIdx-1)*Phi + Qd;
    %% Measurement Update

    % Calculating Range from SV to User based on TRUTH states
    truthPos = [truthStates(1,timeIdx) truthStates(3,timeIdx) truthStates(5,timeIdx)];
    truthVel = [truthStates(2,timeIdx) truthStates(4,timeIdx) truthStates(6,timeIdx)];
    [truthRange,truthRangeRate,~,~] = calcRange(truthPos,truthVel,svPos,svVel);

    % Calculating Estimated Pseudorange from SV to User based on ESTIMATED states
    estiPos = [x_m(1,timeIdx) x_m(3,timeIdx) x_m(5,timeIdx)];
    estiVel = [x_m(2,timeIdx) x_m(4,timeIdx) x_m(6,timeIdx)];
    [range,rangeRate,unitVectors,activeSVs] = calcRange(estiPos,estiVel,svPos,svVel);
    psr = range + x_m(7,timeIdx) - c*satStates(7,:)';
    psrDot = rangeRate - x_m(8,timeIdx);


    j = sum(activeSVs);

    % Generating Correlators, Discriminators, Measurements, and Variances
    [resPsr,resPsrDot,varPsr,varCarr,newCN0,newAmplitude,newNoise,newPhase] = ...
        genCorrelatorResiduals(truthRange(activeSVs),psr(activeSVs),...
        truthRangeRate(activeSVs),psrDot(activeSVs),...
        oldCN0,oldAmplitude,oldNoise,activeSVs,oldPhase,initCN0);

    % Form Z
    count = 1;
    for i = 1:2:2*j
        Z(i:i+1,1) = [resPsr(count);resPsrDot(count)];
        count = count + 1;
    end

    % Form H
    ux = unitVectors(activeSVs,1);
    uy = unitVectors(activeSVs,2);
    uz = unitVectors(activeSVs,3);
    count = 1;
    for i = 1:2:2*j
        H(i:i+1,1:8) = [ux(count) 0 uy(count) 0 uz(count) 0 -1 0;...
            0 ux(count) 0 uy(count) 0 uz(count) 0 -1];
        count = count + 1;
    end

    % Form R
    count = 1;
    Rtmp = [];
    for i = 1:2:2*j
        Rtmp = [Rtmp varPsr(count) varCarr(count)];
        count = count + 1;

    end
    R = diag(Rtmp);

    % Form L
    L = p_m(:,:,timeIdx)*H'*(H*p_m(:,:,timeIdx)*H' + R)^-1;

    % Update States
    x_p = x_m(:,timeIdx) + L*Z;

    % Update Covariance
    p_p = (eye(size(p_m(:,:,timeIdx))) - L*H)*p_m(:,:,timeIdx);

    p_m(:,:,timeIdx) = p_p;
    x_m(:,timeIdx) = x_p;
end

% Position Figures
figure('Position',[200 200 900 800])
hold on
tiledlayout(3,1)
nexttile
hold on
title('ECEF Position Errors')
plot(time,x_m(1,:) - truthStates(1,:),LineWidth=2)
ylabel('ECEF - X [m]')
ax = gca;
ax.FontSize = 16;
nexttile
hold on
plot(time,x_m(3,:) - truthStates(3,:),LineWidth=2)
ylabel('ECEF - Y [m]')
ax = gca;
ax.FontSize = 16;
nexttile
hold on
plot(time,x_m(5,:) - truthStates(5,:),LineWidth=2)
ylabel('ECEF - Z [m]')
xlabel('Time [s]')
ax = gca;
ax.FontSize = 16;


% Velocity Figures
figure('Position',[200 200 900 800])
tiledlayout(3,1)
nexttile
hold on
title('ECEF Velocities')
plot(time,x_m(2,:),LineWidth=2)
plot(time,truthStates(2,:),LineWidth=2)
legend('Estimated','Truth')
ylabel('ECEF - VX [m/s]')
ax = gca;
ax.FontSize = 16;
nexttile
hold on
plot(time,x_m(4,:),LineWidth=2)
plot(time,truthStates(4,:),LineWidth=2)
legend('Estimated','Truth')
ylabel('ECEF - VY [m/s]')
ax = gca;
ax.FontSize = 16;
nexttile
hold on
plot(time,x_m(6,:),LineWidth=2)
plot(time,truthStates(6,:),LineWidth=2)
legend('Estimated','Truth')
ylabel('ECEF - VZ [m/s]')
xlabel('Time [s]')
ax = gca;
ax.FontSize = 16;


% Clock Terms
figure('Position',[200 200 900 800])
hold on
tiledlayout(2,1)
nexttile
hold on
title('Clock Terms')
plot(time,x_m(7,:),LineWidth=2)
plot(time,truthStates(7,:),LineWidth=2)
legend('Estimated','Truth')
ylabel('Clock Bias [m]')
ax = gca;
ax.FontSize = 16;
nexttile
hold on
plot(time,x_m(8,:),LineWidth=2)
plot(time,truthStates(8,:),LineWidth=2)
legend('Estimated','Truth')
ylabel('Clock Drift [m/s]')
ax = gca;
ax.FontSize = 16;

% Geoplot
truthLLA = ecef2lla([truthStates(1,:);truthStates(3,:);truthStates(5,:)]');
estiLLA = ecef2lla([x_m(1,:);x_m(3,:);x_m(5,:)]');

figure('Position',[200 200 900 800])
geoplot(estiLLA(:,1),estiLLA(:,2),'*r')
hold on
geoplot(truthLLA(:,1),truthLLA(:,2),'*k')
legend('Estimated','Truth')
ax = gca;
ax.FontSize = 16;
