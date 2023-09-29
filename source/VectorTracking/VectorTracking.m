function [estimatedStates,refStates,residualPsr,residualCarr,variancePsr,varianceCarr,estimatedCovariance,newCN0,newAmplitude,newNoise,newPhase] = ...
    VectorTracking(forces,moments,predictedStates,receiverStates,predictedCovariance,rinexFilePath,time,year,month,day,oldCN0,oldAmplitude,oldNoise,oldPhase,initCN0,Time_Step,S)

residualPsr = nan(1,31);
residualCarr = nan(1,31);
variancePsr = nan(1,31);
varianceCarr = nan(1,31);
newCN0 = oldCN0';
newAmplitude = oldAmplitude';
newNoise = oldNoise';
newPhase = oldPhase';
estimatedStates = zeros(14,1);
refStates = nan(8,1);

%% Measurement Update
update = 1;
if mod(time,1/50) == 0 && update
    try
        % Convert Receiver States to ECEF
        refStates(1:3) = lla2ecef(receiverStates(7:9)'.*(180/pi),'WGS84') + randn(1,3).*[1.5 1.5 3.0];
        refStates(4:6) = ned2ecefv(receiverStates(1),receiverStates(2),receiverStates(3),receiverStates(7),receiverStates(8),'radians') + randn(1,3).*[0.15 0.15 0.30];
        refStates(7) = 0;
        refStates(8) = 0;

        % Convert Predicted States to ECEF
        estStates(1:3,1) = lla2ecef(predictedStates(7:9)'.*(180/pi),'WGS84');
        estStates(4:6,1) = ned2ecefv(predictedStates(1),predictedStates(2),predictedStates(3),predictedStates(7),predictedStates(8),'radians');
        estStates(7,1) = predictedStates(13);
        estStates(8,1) = predictedStates(14);

        % Pull Current Satellite States
        svStates = genSatellitesStates(time,year,month,day,rinexFilePath);

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
        [resPsr,resCarr,varPsr,varCarr,CN0,Amplitude,Noise,Phase] = genCorrelatorResiduals(refPsr,estPsr,refCarr,estCarr,oldCN0,oldAmplitude,oldNoise,activeSVs,oldPhase,initCN0);

        % Form Z Array
        Z = formZ(resPsr,resCarr);

        % Form H Matrix

        H = formH(unitVectors,predictedStates(7:9));

        % Form R Matrix
        R = formR(varPsr,varCarr);

        % Update Kalman Gain
        L = calcL(H,predictedCovariance,R);

        % Update Estimated States
        updatedStates = updateState(Z,L,predictedStates);
        estimatedStates = updatedStates;

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

    %% Time Update

    % Form Discrete State Transition Matrix (Phi)
    Phi = formPHI(predictedStates,forces,moments,Time_Step);

    % Form Noise Distribution Matrix (Gamma)
    Gamma = formG(Time_Step);

    % Form Disturbance Vector (w)
    w = formW(S,Time_Step);

    % Predict New State
    estimatedStates = Phi*predictedStates + Gamma*w;

    % Form Q
    Qd = formQ(w,Gamma);

    % Predict New Covariance
    estimatedCovariance = Phi*predictedCovariance*Phi' + Qd;

end
end


