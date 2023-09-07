classdef ReceiverPositions < handle
    %RECEIVERPOSITIONS Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = public)
        Property1
    end

    properties (Access = private)
        startTime;
        stopTime;
        timeStep;
        engineParameters;
        WaypointFollower;
        aircraft;

    end

    methods
        function obj = ReceiverPositions(startTime,stopTime,timestep,waypoints,aircraft)

            obj.startTime = startTime;
            obj.stopTime = stopTime;
            obj.timeStep = timestep;
            obj.engineParameters.propR = 400;
            obj.engineParameters.oldFuelFlow = 0;
            obj.engineParameters.oldShaftPower = 0;
            obj.WaypointFollower = uavWaypointFollower("UAVType","fixed-wing","StartFrom","first","Waypoints",waypoints,"TransitionRadius",50);
            obj.aircraft = aircraft;

        end

        function [rcvrStates,controls] = genRCVRStates(obj,initialStates)

            timeArray = obj.startTime:obj.timeStep:obj.stopTime;
            states = initialStates;

            utilities.printText.options(5);
            upd = utilities.progressbar.textprogressbar(length(timeArray),'startmsg','');
            for i = 1:length(timeArray)
                upd(i)

                [augmentedStates,controls] = obj.FVDM_Truth(states);
                states = augmentedStates(4:end);

                rcvrStates(:,i) = obj.body2ecef(augmentedStates);

            end

        end

        function [augmentedState,controls] = FVDM_Truth(obj,state)
            propR = obj.engineParameters.propR;
            oldFuelFlow = obj.engineParameters.oldFuelFlow;
            oldShaftPower = obj.engineParameters.oldShaftPower;

            Vehicle = obj.aircraft.Vehicle;
            BSFC_LUT = obj.aircraft.BSFC_LUT;
            STGeometry = obj.aircraft.STGeometry;
            refLLA = obj.aircraft.LLA;

            % Calculate Atmospheric Parameters for current time step
            LLA = flat2lla(state(7:9)',refLLA(1:2),0,0);
            atmos = Environment('Truth',LLA(1),LLA(2),LLA(3),[],[]);

            % Calculate Controller Reference Values
            [lookaheadPoint,desiredCourse,~,~,~,~] = obj.WaypointFollower([state(7:9);state(12)],obj.aircraft.lookaheadDist);
            desiredCourse = wrapTo2Pi(desiredCourse);
            % Calculate Controller Commands
            % latStick = calcLateralStickCommand(desiredCourse,state);
            % longStick = calcLongitudinalStickCommand(lookaheadPoint,state);
            % pedalPosn = calcPedalPosnCommand(atmos,state);
            % throttle = calcThrottleCommand(70,state);

            latStick = 0;
            longStick = 0;
            pedalPosn = 0;
            throttle = 1;

            controls = [latStick,longStick,pedalPosn,throttle];


            % Propagate Engine and Prop Forces/Moments
            [~,engineForces,engineMoments,propR,oldFuelFlow,oldShaftPower] = PropagateEngine(atmos,state,throttle,propR,obj.timeStep,BSFC_LUT,oldFuelFlow,oldShaftPower);

            % Propagate Aerodynamic Forces/Moments
            [~,aeroForces,aeroMoments] = PropagateAero(atmos.density,state,latStick,longStick,pedalPosn,STGeometry);

            % Propagate Gravity Forces/Moments
            [~,gravityForces,gravityMoments] = PropagateGravity(state,Vehicle.MassProp.Mass,Vehicle.MassProp.r_cg);

            % Summing All Forces and Moments (Body Frame)
            f_ib_b = engineForces + aeroForces + gravityForces;
            m_ib_b = engineMoments + aeroMoments + gravityMoments;

            % Propgating States
            [~,X_dot] = PropagateStates(f_ib_b,m_ib_b,Vehicle.MassProp.MOI,Vehicle.MassProp.r_cg,Vehicle.MassProp.InvMassMat,Vehicle.MassProp.Mass,state,obj.timeStep);
            state = X_dot'*obj.timeStep + state;
            augmentedState = [X_dot(1:3)';state];


            obj.engineParameters.propR = propR;
            obj.engineParameters.oldFuelFlow = oldFuelFlow;
            obj.engineParameters.oldShaftPower = oldShaftPower;
        end

        function rcvrStates = body2ecef(obj,states)

            DCM_b_n = genDCM('rad',[states(15) states(14) states(13)],[3 2 1]);

            pos_NED = DCM_b_n*states(10:12);
            vel_NED = DCM_b_n*states(4:6);
            acc_NED = DCM_b_n*states(1:3);

            [x,y,z] = ned2ecef(pos_NED(1),pos_NED(2),pos_NED(3),obj.aircraft.LLA(1),obj.aircraft.LLA(2),0,wgs84Ellipsoid("meter"),"degrees");
            [u,v,w] = ned2ecefv(vel_NED(1),vel_NED(2),vel_NED(3),obj.aircraft.LLA(1),obj.aircraft.LLA(2),"degrees");
            [ax,ay,az] = ned2ecefv(acc_NED(1),acc_NED(2),acc_NED(3),obj.aircraft.LLA(1),obj.aircraft.LLA(2),"degrees");

            rcvrStates = [x,u,ax,y,v,ay,z,w,az,0,0];

        end
    end
end

