%% Formatting
clc
clear
close all
format shortg

%% Static Vector Tracking Analysis
AuburnLLA = [32.6099,85.4808,200];
AuburnECEF = lla2ecef(AuburnLLA,'WGS84');
x_m(:,1) = [AuburnECEF(1);0;AuburnECEF(2);0;AuburnECEF(3);0;0;0];
timeStep = 1/50;
k = [0 1;0 0];
A = blkdiag(k,k,k,k);
Phi = expm(A*timeStep);




for i = 2:100

    satStates = 1;
    j = length(satStates(1,:));

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