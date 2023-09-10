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
        function obj = ReceiverPositions(startTime,stopTime,timestep,aircraft)

            obj.startTime = startTime;
            obj.stopTime = stopTime;
            obj.timeStep = timestep;
            obj.engineParameters.propR = 400;
            obj.engineParameters.oldFuelFlow = 0;
            obj.engineParameters.oldShaftPower = 0;
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

        function state = noiseFVDM(obj,state,controls)
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

            throttle = controls(4);
            longStick = controls(2);
            latStick = controls(1);
            pedalPosn = controls(3);

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
            state = X_dot*obj.timeStep + state;


            obj.engineParameters.propR = propR;
            obj.engineParameters.oldFuelFlow = oldFuelFlow;
            obj.engineParameters.oldShaftPower = oldShaftPower;
        end

        function rcvrStates = body2ecef(obj,states)

            DCM_b_n = genDCM('rad',[states(12) states(11) states(10)],[3 2 1]);

            pos_NED = DCM_b_n*states(7:9);
            vel_NED = DCM_b_n*states(1:3);
        

            [x,y,z] = ned2ecef(pos_NED(1),pos_NED(2),pos_NED(3),obj.aircraft.LLA(1),obj.aircraft.LLA(2),0,wgs84Ellipsoid("meter"),"degrees");
            [u,v,w] = ned2ecefv(vel_NED(1),vel_NED(2),vel_NED(3),obj.aircraft.LLA(1),obj.aircraft.LLA(2),"degrees");
     

            rcvrStates = [x,u,y,v,z,w,0,0];

        end

        function cov = covariance(obj,p_m,bodyStates)
            stop = 1;

        end
    end
end

