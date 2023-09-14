function [estimatedStates,refStates,estimatedCovariance] = VectorTracking(predictedStates,receiverStates,predictedCovariance,satelliteStates,time,refLLA)

Q  = diag([0.1 0.1 0.1 1e-6 1e-6 1e-6 1 1 1 1e-5 1e-5 1e-5 1e-7 1e-8]);
timeStep = 1/400;

%% Time Update
estimatedStates = predictedStates;

PHI = calcJacobian(predictedStates,timeStep);

estimatedCovariance = PHI*predictedCovariance*PHI' + Q;

refStates = receiverStates;

%% Measurement Update

if mod(time,1/50) == 0

    % Convert Receiver States to ECEF
    refStates = flat2ecef(refLLA,receiverStates);
    refStates = refStates + randn(12,1).*[0.15;0.15;0.15;0;0;0;1.5;1.5;3.0;0;0;0];

    % Convert Predicted States to ECEF
    estStates = flat2ecef(refLLA,predictedStates);

    % Pull Current Satellite States
    svStates = satelliteStates(:,:,floor(time*50 + 1));

    % Calculate Reciever PSR, Carrier Frequency, Unit Vectors
    [refPsr,refCarr,unitVectors] = calcPsr(refStates,svStates);

    % Calculate Predicted PSR, Carrier Frequency
    [estPsr,estCarr,~] = calcPsr(estStates,svStates);

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
end

end

