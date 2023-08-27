function [state,engineParameters,controls] = FVDM_Truth(deltaTime,state,engineParameters,WaypointFollower,aircraft)
propR = engineParameters.propR;
oldFuelFlow = engineParameters.oldFuelFlow;
oldShaftPower = engineParameters.oldShaftPower;

Vehicle = aircraft.Vehicle;
BSFC_LUT = aircraft.BSFC_LUT;
STGeometry = aircraft.STGeometry;
refLLA = aircraft.LLA;

% Calculate Atmospheric Parameters for current time step
LLA = flat2lla(state(7:9)',refLLA(1:2),0,0);
atmos = Environment('Truth',LLA(1),LLA(2),LLA(3),[],[]);

% Calculate Controller Reference Values
[lookaheadPoint,desiredCourse,~,~,~,~] = WaypointFollower([state(7:9);state(12)],aircraft.lookaheadDist);
desiredCourse = wrapTo2Pi(desiredCourse);
% Calculate Controller Commands
latStick = calcLateralStickCommand(desiredCourse,state);
longStick = calcLongitudinalStickCommand(lookaheadPoint,state);
pedalPosn = calcPedalPosnCommand(atmos,state);
throttle = calcThrottleCommand(70,state);

controls = [latStick,longStick,pedalPosn,throttle];


% Propagate Engine and Prop Forces/Moments
[~,engineForces,engineMoments,propR,oldFuelFlow,oldShaftPower] = PropagateEngine(atmos,state,throttle,propR,deltaTime,BSFC_LUT,oldFuelFlow,oldShaftPower);

% Propagate Aerodynamic Forces/Moments
[~,aeroForces,aeroMoments] = PropagateAero(atmos.density,state,latStick,longStick,pedalPosn,STGeometry);

% Propagate Gravity Forces/Moments
[~,gravityForces,gravityMoments] = PropagateGravity(state,Vehicle.MassProp.Mass,Vehicle.MassProp.r_cg);

% Summing All Forces and Moments (Body Frame)
f_ib_b = engineForces + aeroForces + gravityForces;
m_ib_b = engineMoments + aeroMoments + gravityMoments;

% Propgating States
[~,state] = PropagateStates(f_ib_b,m_ib_b,Vehicle.MassProp.MOI,Vehicle.MassProp.r_cg,Vehicle.MassProp.InvMassMat,Vehicle.MassProp.Mass,state,deltaTime);

engineParameters.propR = propR;
engineParameters.oldFuelFlow = oldFuelFlow;
engineParameters.oldShaftPower = oldShaftPower;
end