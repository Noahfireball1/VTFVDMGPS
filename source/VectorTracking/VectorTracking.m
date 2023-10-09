function [estimatedStates,receiverStates,residualPsr,residualCarr,variancePsr,varianceCarr,estimatedCovariance,newCN0,newAmplitude,newNoise,newPhase,satStates] = ...
    VectorTracking(forces,moments,predictedStates,trueF_ib_b,trueM_ib_b,receiverStates,predictedCovariance,rinexFilePath,time,year,month,day,oldCN0,oldAmplitude,oldNoise,...
    oldPhase,initCN0,Time_Step,variance,clkVar,Q)

residualPsr = nan(1,31);
residualCarr = nan(1,31);
variancePsr = nan(1,31);
varianceCarr = nan(1,31);
newCN0 = oldCN0';
newAmplitude = oldAmplitude';
newNoise = oldNoise';
newPhase = oldPhase';
satStates = zeros(7,31);

c = 299792458;
e = 0.0818191910428; % eccentricity
R_0 = 6378137.0; % equatorial radius [meters]

%% Propagate Receiver States
truthStates_n = equationsOfMotion(receiverStates,trueF_ib_b,trueM_ib_b,Time_Step);
truthStates_n = truthStates_n; %+ sqrt(Q)*randn(14,1)*Time_Step;
% For Conversion to ECEF Frame
trueLat = truthStates_n(7);
trueLong = truthStates_n(8);
trueAlt = truthStates_n(9);

% Matrix Rotating from Navigation to ECEF Frame (Groves Eq. 2.150)
trueC_n_e = [-cos(trueLong)*sin(trueLat) -sin(trueLong)*sin(trueLat) cos(trueLat);...
    -sin(trueLong)                     cos(trueLong)           0      ;...
    -cos(trueLong)*cos(trueLat)  -cos(trueLat)*sin(trueLong) -sin(trueLat)]';

% Tranverse Radius of Curvature (Groves Eq. 2.106)
trueR_E = (R_0)/(sqrt(1 - (e^2)*sin(trueLat)^2));

%% Time Update

% Predict New State
estStates_n = equationsOfMotion(predictedStates,forces,moments,Time_Step);

% For Conversion to ECEF Frame
estiLat = estStates_n(7);
estiLong = estStates_n(8);
estiAlt = estStates_n(9);

% Matrix Rotating from Navigation to ECEF Frame (Groves Eq. 2.150)
estiC_n_e = [-cos(estiLong)*sin(estiLat) -sin(estiLong)*sin(estiLat) cos(estiLat);...
    -sin(estiLong)                     cos(estiLong)           0      ;...
    -cos(estiLong)*cos(estiLat)  -cos(estiLat)*sin(estiLong) -sin(estiLat)]';

% Tranverse Radius of Curvature (Groves Eq. 2.106)
estiR_E = (R_0)/(sqrt(1 - (e^2)*sin(estiLat)^2));

% Form Discrete State Transition Matrix (Phi)
Phi = formPHI(predictedStates,forces,moments,Time_Step);

% Predict New Covariance
predictedCovariance = Phi*predictedCovariance*Phi' + Phi*blkdiag(Q(1:12,1:12)*Time_Step,clkVar)*Phi';

%% Measurement Update
update = 1;
% if time > 60
%     initCN0 = 10^(22.5/10)*ones(1,31);
% end
if mod(time,1/50) == 0 && update
    try
        %% Convert True States to ECEF Frame [Positions; Velocities; Clock Bias; Clock Drift]

        % Position Groves (Eq. 2.112)
        truthStates_e(1:3) = [(trueR_E + trueAlt)*cos(trueLat)*cos(trueLong);...
            (trueR_E + trueAlt)*cos(trueLat)*sin(trueLong);...
            ((1 - e^2)*trueR_E + trueAlt)*sin(trueLat)];

        % Velocity Groves (Eq. 2.150)
        truthStates_e(4:6) = trueC_n_e*truthStates_n(1:3);

        % Clock Terms
        truthStates_e(7:8) = truthStates_n(13:14);

        %% Convert Estimated States to ECEF Frame [Positions; Velocities; Clock Bias; Clock Drift]

        % Position Groves (Eq. 2.112)
        estStates_e(1:3) = [(estiR_E + estiAlt)*cos(estiLat)*cos(estiLong);...
            (estiR_E + estiAlt)*cos(estiLat)*sin(estiLong);...
            ((1 - e^2)*estiR_E + estiAlt)*sin(estiLat)];

        % Velocity Groves (Eq. 2.150)
        estStates_e(4:6) = estiC_n_e*estStates_n(1:3);

        % Clock Terms
        estStates_e(7:8) = estStates_n(13:14);

        %% Pull Current Satellite States
        satStates = genSatellitesStates(time,year,month,day,rinexFilePath);
        svPos = [satStates(1,:)' satStates(3,:)' satStates(5,:)'];
        svVel = [satStates(2,:)' satStates(4,:)' satStates(6,:)'];

        % Calculating Range from SV to User based on TRUTH states
        truthPos = [truthStates_e(1) truthStates_e(2) truthStates_e(3)];
        truthVel = [truthStates_e(4) truthStates_e(5) truthStates_e(6)];
        [truthRange,truthRangeRate,~,~] = calcRange(truthPos,truthVel,svPos,svVel);

        % Calculating Estimated Pseudorange from SV to User based on ESTIMATED states
        estiPos = [estStates_e(1) estStates_e(2) estStates_e(3)];
        estiVel = [estStates_e(4) estStates_e(5) estStates_e(6)];
        [range,rangeRate,unitVectors,activeSVs] = calcRange(estiPos,estiVel,svPos,svVel);
        psr = range + estStates_e(7) - c*satStates(7,:)';
        psrDot = rangeRate - estStates_e(8);

        % Generate Correlator Residuals
        [resPsr,resCarr,varPsr,varCarr,CN0,Amplitude,Noise,Phase] = genCorrelatorResiduals(truthRange(activeSVs),psr(activeSVs),truthRangeRate(activeSVs),psrDot(activeSVs),...
            oldCN0,oldAmplitude,oldNoise,activeSVs,oldPhase,initCN0);

        % Form Z Array
        Z = formZ(resPsr,resCarr);

        % Form H Matrix
        H = formH(unitVectors(activeSVs,:),estStates_n(7:9));

        % Form R Matrix
        R = formR(varPsr,varCarr);

        % Form L
        L = predictedCovariance*H'*(H*predictedCovariance*H' + R)^-1;

        % Update Estimated States
        estimatedStates = estStates_n + L*Z;
        % estimatedStates = estStates_n;

        % Update Estimated Covariance
        estimatedCovariance = (eye(size(predictedCovariance)) - L*H)*predictedCovariance;

        % Convert residuals to determinstic size
        activeIdx = find(activeSVs);
        residualPsr(activeIdx) = resPsr;
        residualCarr(activeIdx) = resCarr;
        variancePsr(activeIdx) = varPsr;
        varianceCarr(activeIdx) = varCarr;
        newCN0(activeIdx) = CN0;
        newAmplitude(activeIdx) = Amplitude;
        newNoise(activeIdx) = Noise;
        newPhase(activeIdx) = Phase;
    catch
        sprintf('%f',time)

    end

else
    estimatedStates = estStates_n;
    estimatedCovariance = predictedCovariance;

end
receiverStates = truthStates_n;

end

