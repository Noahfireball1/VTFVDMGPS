function [estimatedStates,refStates,residualPsr,residualCarr,variancePsr,varianceCarr,estimatedCovariance] = VectorTracking(predictedStates,receiverStates,predictedCovariance,rinexFilePath,time,year,month,day,refLLA)

Q  = diag([0.1 0.1 0.1 1e-6 1e-6 1e-6 10 10 10 1e-5 1e-5 1e-5 1e-7 1e-8]);
% Q  = zeros(14);
timeStep = 1/400;

%% Time Update
estimatedStates = predictedStates;

PHI = calcJacobian(predictedStates,timeStep);

estimatedCovariance = PHI*predictedCovariance*PHI' + Q;

refStates = nan(12,1);
residualPsr = nan(1,31);
residualCarr = nan(1,31);
variancePsr = nan(1,31);
varianceCarr = nan(1,31);

%% Measurement Update

if mod(time,1/50) == 0

    % Convert Receiver States to ECEF
    refStates = flat2ecef(refLLA,receiverStates);
    refStates = refStates + randn(12,1).*[0.15;0.15;0.15;0;0;0;1.5;1.5;3.0;0;0;0];

    % Convert Predicted States to ECEF
    estStates = flat2ecef(refLLA,predictedStates);

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
    [resPsr,resCarr,varPsr,varCarr] = genCorrelatorResiduals(refPsr,estPsr,refCarr,estCarr);

    % Form Z Array
    Z = formZ(resPsr,resCarr);

    % Form H Matrix
    H = formH(unitVectors);

    % Form R Matrix
    R = formR(varPsr,varCarr);

    % Update Kalman Gain
    L = calcL(H,predictedCovariance,R);

    % Update Estimated States
    estStates = updateState(Z,L,estStates);

    % Update Estimated Covariance
    estimatedCovariance = updateCovariance(predictedCovariance,L,H);

    % Convert States Back to Flat Earth Frame
    estimatedStates = ecef2flat(refLLA,estStates);
    refStates = ecef2flat(refLLA,refStates);

    % Convert residuals to determinstic size
    activeIdx = find(activeSVs);
    residualPsr(activeIdx) = resPsr; 
    residualCarr(activeIdx) = resCarr;
    variancePsr(activeIdx) = varPsr;
    varianceCarr(activeIdx) = varCarr;
end

end

