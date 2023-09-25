function [estimatedStates,refStates,residualPsr,residualCarr,variancePsr,varianceCarr,estimatedCovariance,newCN0,newAmplitude,newNoise,newPhase] = ...
    VectorTracking(predictedStates,receiverStates,predictedCovariance,rinexFilePath,time,year,month,day,refLLA,timeStep,oldCN0,oldAmplitude,oldNoise,m,cg,imm,oldPhase,initCN0)

residualPsr = nan(1,31);
residualCarr = nan(1,31);
variancePsr = nan(1,31);
varianceCarr = nan(1,31);
newCN0 = oldCN0';
newAmplitude = oldAmplitude';
newNoise = oldNoise';
newPhase = oldPhase';
estimatedStates = predictedStates;
refStates = nan(12,1);

%% Measurement Update
if mod(time,1/50) == 0
    try
        % Convert Receiver States to ECEF
        refStates = receiverStates;

        % + randn(12,1).*[0.15;0.15;0.15;0;0;0;1.5;1.5;3.0;0;0;0];

        % Convert Predicted States to ECEF
        estStates = predictedStates;
        
        % Pull Current Satellite States
        svStates = genSatellitesStates(time,year,month,day,rinexFilePath);

        % Calculate Reciever PSR, Carrier Frequency, Unit Vectors
        [psr,carr,unitVectors,activeSVs] = calcPsr(refStates,svStates);
        refPsr = psr(activeSVs);
        refCarr = carr(activeSVs);
        unitVectors = unitVectors(:,activeSVs);

        % Calculate Predicted PSR, Carrier Frequency
        [psr,carr,~,~] = calcPsr(estStates,svStates);
        estPsr = psr(activeSVs);
        estCarr = carr(activeSVs);

        % Generate Correlator Residuals
        [resPsr,resCarr,varPsr,varCarr,CN0,Amplitude,Noise,Phase] = genCorrelatorResiduals(refPsr,estPsr,refCarr,estCarr,oldCN0,oldAmplitude,oldNoise,activeSVs,oldPhase,initCN0);

        % Form Z Array
        Z = formZ(resPsr,resCarr);

        % Form H Matrix

        H = formH(unitVectors);

        % Form R Matrix
        R = formR(varPsr,varCarr);

        % Update Kalman Gain
        L = calcL(H,predictedCovariance,R);

        % Update Estimated States
        updStates = updateState(Z,L,estStates);

        % Update Estimated Covariance
        estimatedCovariance = updateCovariance(predictedCovariance,L,H);

        % Convert States Back to Flat Earth Frame
        estimatedStates = ecef2flat(refLLA,updStates);
        refStates = ecef2flat(refLLA,refStates);

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

    estimatedStates = predictedStates;

    % Update Covariance
    [estimatedCovariance] = calculateCovariance(predictedStates,timeStep,predictedCovariance,refLLA);

end
end


