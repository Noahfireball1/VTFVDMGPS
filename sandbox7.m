%% Formatting
clc
clear
close all
format shortg

rinexFilePath = "C:\Users\nsm0014\graduate\thesis\VTFVDMGPS\source\VectorTracking\Satellites\rinex\AMC400USA_R_20223560000_01D_GN.rnx";
%% Static Vector Tracking Analysis
AuburnLLA = [32.6099,85.4808,200];
AuburnECEF = lla2ecef(AuburnLLA,'WGS84');
x_m(:,1) = [AuburnECEF(1);0;AuburnECEF(2);0;AuburnECEF(3);0;0;0];
timeStep = 1/50;
time = 0:timeStep:100;
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


for i = 2:100

    satStates = genSatellitesStates(time(i),2022,12,22,rinexFilePath);
    %% Truth State Propagation
    truthStates(:,i) = Phi*truthStates(:,i-1) + sqrt(Qd)*randn(8,1);

    %% Time Update
    x_m(:,i) = Phi*x_m(:,i-1);
    p_m(:,:,i) = Phi*p_m(:,:,i-1)*Phi + Qd;


    %% Measurement Update
    trueRange = 1;
    trueRangeRate = 1;

    estiRange = 1;
    estiRangeRate = 1;

    % Calculate Code Phase Error
    codePhaseError = 1;

    % Calculate Carrier Frequency Error
    carrierFreqError = 1;

    % Generate Correlators
    for i = 1:length(range)

    end

    % Caclulate Discriminators
    dllDisc = 1;
    fllDisc = 1;

    % Calculate Noise Terms
    gpsNoise = 1;
    newAmplitude = 1;
    newNoise = 1;

    % Calculate Measurement Variance
    varPsr = 1;
    varCarr = 1;

    % Form Z
    z = [resPsr;resPsrDot];
    Z = repmat(z,j,1);

    % Form H
    h = [ux 0 uy 0 uz 0 1 0;...
         0 ux 0 uy 0 uz 0 1];
    H = repmat(h,j,1);

    % Form R
    r = diag([varPsr,varCarr]);
    R = repmat(r,j,1);

    % Form L
    L = p_m*H'*(H*p_m*H' + R)^-1;

    % Update States
    x_p = x_m - L*Z;

    % Update Covariance
    p_p = (eye(size(p_m)) - L*H)*p_m;
end