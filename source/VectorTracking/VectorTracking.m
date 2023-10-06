function [estimatedStates,receiverStates,refStates,residualPsr,residualCarr,variancePsr,varianceCarr,estimatedCovariance,newCN0,newAmplitude,newNoise,newPhase,svStates] = ...
    VectorTracking(forces,moments,predictedStates,trueF_ib_b,trueM_ib_b,receiverStates,predictedCovariance,rinexFilePath,time,year,month,day,oldCN0,oldAmplitude,oldNoise,oldPhase,initCN0,Time_Step,variance,clkVar)

residualPsr = nan(1,31);
residualCarr = nan(1,31);
variancePsr = nan(1,31);
varianceCarr = nan(1,31);
newCN0 = oldCN0';
newAmplitude = oldAmplitude';
newNoise = oldNoise';
newPhase = oldPhase';
refStates = nan(8,1);
svStates = zeros(7,31);

%% Propagate Receiver States
receiverStates = truthModelStateUpdate(trueF_ib_b,trueM_ib_b,receiverStates,Time_Step,variance,clkVar);

%% Time Update
% Form Discrete State Transition Matrix (Phi)
% Calculate Q
Q = blkdiag(diag([variance(1),variance(2),variance(3)...
    ,0.001,0.001,0.001 ...
    ,0,0,0 ...
    ,0.001 0.001 0.001])...
    ,clkVar);
% Q = blkdiag(diag([variance(1) variance(2), variance(3)]),zeros(9),clkVar);

% Predict New State
predictedStates = predictStates(predictedStates,forces,moments,Time_Step);
Phi = formPHI(predictedStates,forces,moments,Time_Step);
% Predict New Covariance
predictedCovariance = Phi*predictedCovariance*Phi' + Phi*Q*Phi';

%% Measurement Update
update = 1;
if mod(time,1/50) == 0 && update
    try
        % Convert Receiver States to ECEF

        refStates(1:3) = lla2ecef(receiverStates(7:9)'.*(180/pi),'WGS84')';
        [refStates(4,1),refStates(5,1),refStates(6,1)] = ned2ecefv(receiverStates(1),receiverStates(2),receiverStates(3),receiverStates(7),receiverStates(8),'radians');
        refStates(7) = receiverStates(13);
        refStates(8) = receiverStates(14);

        % Convert Predicted States to ECEF
        estStates(1:3,1) = lla2ecef(predictedStates(7:9)'.*(180/pi),'WGS84');
        [estStates(4,1),estStates(5,1),estStates(6,1)] = ned2ecefv(predictedStates(1),predictedStates(2),predictedStates(3),predictedStates(7),predictedStates(8),'radians');
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

        % Generate Correlator Residuals
        [resPsr,resCarr,varPsr,varCarr,CN0,Amplitude,Noise,Phase] = genCorrelatorResiduals(refPsr,estPsr,refCarr,estCarr,oldCN0,oldAmplitude,oldNoise,activeSVs,oldPhase,initCN0);

        % Form Z Array
        Z = formZ(resPsr,resCarr);

        % Form H Matrix
        H = formH(unitVectors(:,activeSVs),predictedStates(7:9));

        % Form R Matrix
        R = formR(varPsr,varCarr);

        % Update Kalman Gain
        L = calcL(H,predictedCovariance,R);

        % Update Estimated States
        estimatedStates = predictedStates - L*Z;

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
        if time == 13
            stop = 1;
        end
    catch
        sprintf('%f',time)

    end

else
    estimatedStates = predictedStates;
    estimatedCovariance = predictedCovariance;

end



end

