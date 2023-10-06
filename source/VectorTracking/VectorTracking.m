function [estimatedStates,receiverStates,residualPsr,residualCarr,variancePsr,varianceCarr,estimatedCovariance,newCN0,newAmplitude,newNoise,newPhase,svStates] = ...
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
svStates = zeros(7,31);

e = 0.0818191910428; % eccentricity
R_0 = 6378137.0; % equatorial radius [meters]

%% Propagate Receiver States
refStates_n = equationsOfMotionWithNoise(receiverStates,trueF_ib_b,trueM_ib_b,Time_Step,variance,clkVar);

% For Conversion to ECEF Frame
trueLat = refStates_n(7);
trueLong = refStates_n(8);
trueAlt = refStates_n(9);

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
Phi = formPHI(estStates_n,forces,moments,Time_Step);

% Predict New Covariance
predictedCovariance = Phi*predictedCovariance*Phi' + blkdiag(Q(1:12,1:12)*Time_Step,clkVar);

%% Measurement Update
update = 1;
if mod(time,1/50) == 0 && update
    try
        %% Convert True States to ECEF Frame [Positions; Velocities; Clock Bias; Clock Drift]

        % Position Groves (Eq. 2.112)
        refStates_e(1:3) = [(trueR_E + trueAlt)*cos(trueLat)*cos(trueLong);...
            (trueR_E + trueAlt)*cos(trueLat)*sin(trueLong);...
            ((1 - e^2)*trueR_E + trueAlt)*sin(trueLat)];

        % Velocity Groves (Eq. 2.150)
        refStates_e(4:6) = trueC_n_e*refStates_n(1:3);

        % Clock Terms
        refStates_e(7:8) = refStates_n(13:14);

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
        svStates = genSatellitesStates(time,year,month,day,rinexFilePath);

        % Calculate Reciever PSR, Carrier Frequency, Unit Vectors
        [psr,carr,~,activeSVs] = calcPsr(refStates_e,svStates);
        refPsr = psr(activeSVs);
        refCarr = carr(activeSVs);

        % Calculate Predicted PSR, Carrier Frequency
        [psr,carr,unitVectors,~] = calcPsr(estStates_e,svStates);
        estPsr = psr(activeSVs);
        estCarr = carr(activeSVs);

        % Generate Correlator Residuals
        [resPsr,resCarr,varPsr,varCarr,CN0,Amplitude,Noise,Phase] = genCorrelatorResiduals(refPsr,estPsr,refCarr,estCarr,oldCN0,oldAmplitude,oldNoise,activeSVs,oldPhase,initCN0);

        % Form Z Array
        Z = formZ(resPsr,resCarr);

        % Form H Matrix
        H = formH(unitVectors(:,activeSVs),estStates_n(7:9));

        % Form R Matrix
        R = formR(varPsr,varCarr);

        % Update Kalman Gain
        L = calcL(H,predictedCovariance,R);

        % Update Estimated States
        estimatedStates = estStates_n - L*Z;

        % Update Estimated Covariance
        estimatedCovariance = updateCovariance(predictedCovariance,L,H,R);

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
receiverStates = refStates_n;

end

