function [state,covariance,engineParameters] = FVDM_Noise(deltaTime,state,engineParameters,controls,aircraft)
propR = engineParameters.propR;
oldFuelFlow = engineParameters.oldFuelFlow;
oldShaftPower = engineParameters.oldShaftPower;

Vehicle = aircraft.Vehicle;
BSFC_LUT = aircraft.BSFC_LUT;
STGeometry = aircraft.STGeometry;
refLLA = aircraft.LLA;
engineForcesVAR = aircraft.engineForcesVAR;
engineMomentsVAR = aircraft.engineMomentsVAR;
aeroForcesVAR = aircraft.aeroForcesVAR;
aeroMomentsVAR = aircraft.aeroMomentsVAR;
gravityForcesVAR = aircraft.gravityForcesVAR;
gravityMomentsVAR = aircraft.gravityMomentsVAR;

% Calculate Atmospheric Parameters for current time step
LLA = flat2lla(state(7:9)',refLLA(1:2),0,0);
atmos = Environment('Noise',LLA(1),LLA(2),LLA(3),50,50);

% Calculate Controller Commands
latStick = controls(1);
longStick = controls(2);
pedalPosn = controls(3);
throttle = controls(4);


% Propagate Engine and Prop Forces/Moments
[~,engineForces,engineMoments,propR,oldFuelFlow,oldShaftPower] = PropagateEngine(atmos,state,throttle,propR,deltaTime,BSFC_LUT,oldFuelFlow,oldShaftPower);
engineForces = engineForces + randn(1,3)*engineForcesVAR;
engineMoments = engineMoments + randn(1,3)*engineMomentsVAR;
% engineForces = 0;
% engineMoments = 0;

% Propagate Aerodynamic Forces/Moments
[~,aeroForces,aeroMoments] = PropagateAero(atmos.density,state,latStick,longStick,pedalPosn,STGeometry);
aeroForces = aeroForces + randn(1,3)*aeroForcesVAR;
aeroMoments = aeroMoments + randn(1,3)*aeroMomentsVAR;

% Propagate Gravity Forces/Moments
[~,gravityForces,gravityMoments] = PropagateGravity(state,Vehicle.MassProp.Mass,Vehicle.MassProp.r_cg);
gravityForces = gravityForces + randn(1,3)*gravityForcesVAR;
gravityMoments = gravityMoments + randn(1,3)*gravityMomentsVAR;

% Summing All Forces and Moments (Body Frame)
f_ib_b = engineForces + aeroForces + gravityForces;
m_ib_b = engineMoments + aeroMoments + gravityMoments;

% Propgating States
[~,state,covariance] = PropagateStates(f_ib_b,m_ib_b,Vehicle.MassProp.MOI,Vehicle.MassProp.r_cg,Vehicle.MassProp.InvMassMat,Vehicle.MassProp.Mass,state,deltaTime);

engineParameters.propR = propR;
engineParameters.oldFuelFlow = oldFuelFlow;
engineParameters.oldShaftPower = oldShaftPower;
end