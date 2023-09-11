classdef ReceiverPositions < handle
    %RECEIVERPOSITIONS Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = public)
        Q;
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
            obj.Q = diag([0.01 0.01 0.01 0.0001 0.0001 0.0001 0.1 0.1 0.1 0.001 0.001 0.001 0.0000001 0.000001]);

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

            [x,y,z] = ned2ecef(states(7),states(8),states(9),obj.aircraft.LLA(1),obj.aircraft.LLA(2),0,wgs84Ellipsoid("meter"),"degrees");
            [u,v,w] = ned2ecefv(states(1),states(2),states(3),obj.aircraft.LLA(1),obj.aircraft.LLA(2),"degrees");

            rcvrStates = [u,v,w,states(4:6)',x,y,z,states(10:12)',0,0];

        end

        function NEDStates = ecef2body(obj,ecefState)

            [N,E,D] = ecef2ned(ecefState(7),ecefState(8),ecefState(9),obj.aircraft.LLA(1),obj.aircraft.LLA(2),0,wgs84Ellipsoid('meter'));
            [Nv,Ev,Dv] = ecef2nedv(ecefState(1),ecefState(2),ecefState(3),obj.aircraft.LLA(1),obj.aircraft.LLA(2));

            NEDStates = [Nv,Ev,Dv,ecefState(4:6)',N,E,D,ecefState(10:12)']';

        end

        function cov = covariance(obj,timeStep,state,p_m)

            u = state(1);
            v = state(2);
            w = state(3);
            p = state(4);
            q = state(5);
            r = state(6);
            phi = state(10);
            theta = state(11);
            ps = state(12);

            f1 = [-(298977403082689*q)/1267650600228229401496703205376,...
                (298977403082689*p)/1267650600228229401496703205376 + (7205759403792793197*r)/7205759403792793600,...
                -(7205759403792793197*q)/7205759403792793600,...
                (1067922029527797*p)/45035996273704960 - (13143461160180030981*r)/36028797018963968000 + (298977403082689*v)/1267650600228229401496703205376,...
                - (479633360314957851*q)/562949953421312000 - (298977403082689*u)/1267650600228229401496703205376 - (7205759403792793197*w)/7205759403792793600,...
                (7205759403792793197*v)/7205759403792793600 - (985964771368110627*r)/1125899906842624000 - (13143461160180030981*p)/36028797018963968000,...
                0,...
                0,...
                0,...
                0,...
                0,...
                0];

            f2 = [    -(507060240091291686981746483975017*r)/507060240091291760598681282150400,...
                0,...
                (507060240091291686981746483975017*p)/507060240091291760598681282150400,...
                (176340365164395058951*q)/288230376151711744000 + (507060240091291686981746483975017*w)/507060240091291760598681282150400,...
                (176340365164395058951*p)/288230376151711744000 - (590980669729793734301*r)/1475739525896764129280,...
                - (590980669729793734301*q)/1475739525896764129280 - (507060240091291686981746483975017*u)/507060240091291760598681282150400,...
                0,...
                0,...
                0,...
                0,...
                0,...
                0];

            f3 = [(162259276829213384308593495216085*q)/162259276829213363391578010288128,...
                (963*r)/184467440737095516160 - (162259276829213384308593495216085*p)/162259276829213363391578010288128,...
                -(963*q)/184467440737095516160, (507226810809405*p)/1125899906842624 + (13828929479546163813*r)/184467440737095516160 - (162259276829213384308593495216085*v)/162259276829213363391578010288128,...
                (9007199254740993*q)/22517998136852480 + (162259276829213384308593495216085*u)/162259276829213363391578010288128 - (963*w)/184467440737095516160,...
                (13828929479546163813*p)/184467440737095516160 - (1137336961447107*r)/22517998136852480 + (963*v)/184467440737095516160,...
                0,...
                0,...
                0,...
                0,...
                0,...
                0];

            f4 = [(5824665028096681*r)/3245185536584267267831560205762560,...
                0,...
                -(5824665028096681*p)/3245185536584267267831560205762560,...
                (335585127691199286301*q)/9223372036854775808000 - (5824665028096681*w)/3245185536584267267831560205762560,...
                (335585127691199286301*p)/9223372036854775808000 - (42649887851095123581*r)/46116860184273879040,...
                (5824665028096681*u)/3245185536584267267831560205762560 - (42649887851095123581*q)/46116860184273879040,...
                0,...
                0,...
                0,...
                0,...
                0,...
                0];

            f5 = [(545508014063655*q)/40564819207303340847894502572032,...
                (1253*r)/230584300921369395200 - (545508014063655*p)/40564819207303340847894502572032,...
                -(1253*q)/230584300921369395200,...
                (950032281021957115519*r)/1152921504606846976000 - (1334902536909751*p)/11258999068426240 - (545508014063655*v)/40564819207303340847894502572032,...
                q/140737488355328000 + (545508014063655*u)/40564819207303340847894502572032 - (1253*w)/230584300921369395200,...
                (950032281021957115519*p)/1152921504606846976000 + (33372563422743777*r)/281474976710656000 + (1253*v)/230584300921369395200,...
                0,...
                0,...
                0,...
                0,...
                0,...
                0];

            f6 = [(1911107816070405*r)/324518553658426726783156020576256,...
                0,...
                -(1911107816070405*p)/324518553658426726783156020576256,...
                (1932650436645849243469*q)/4611686018427387904000 - (1911107816070405*w)/324518553658426726783156020576256,...
                (1932650436645849243469*p)/4611686018427387904000 - (429548963444735107403*r)/11805916207174113034240,...
                (1911107816070405*u)/324518553658426726783156020576256 - (429548963444735107403*q)/11805916207174113034240,...
                0,...
                0,...
                0,...
                0,...
                0,...
                0];

            f7 = [ cos(ps)*cos(theta),...
                -cos(phi)*sin(ps),...
                sin(phi)*sin(ps),...
                0,...
                0,...
                0,...
                0,...
                0,...
                0,...
                w*cos(phi)*sin(ps) + v*sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta) - cos(ps)*sin(phi)*sin(theta),...
                cos(phi)*cos(ps)*cos(theta) - u*cos(ps)*sin(theta) + cos(ps)*cos(theta)*sin(phi),...
                w*cos(ps)*sin(phi) - v*cos(phi)*cos(ps) - sin(phi)*sin(ps)*sin(theta) - u*cos(theta)*sin(ps) - cos(phi)*sin(ps)*sin(theta)];

            f8 = [    cos(theta)*sin(ps),...
                cos(phi)*cos(ps),...
                -cos(ps)*sin(phi),...
                0,...
                0,...
                0,...
                0,...
                0,...
                0,...
                cos(phi)*sin(ps)*sin(theta) - w*cos(phi)*cos(ps) - v*cos(ps)*sin(phi) - sin(phi)*sin(ps)*sin(theta),...
                cos(phi)*cos(theta)*sin(ps) - u*sin(ps)*sin(theta) + cos(theta)*sin(phi)*sin(ps),...
                u*cos(ps)*cos(theta) - v*cos(phi)*sin(ps) + w*sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta) + cos(ps)*sin(phi)*sin(theta)];

            f9 = [ -sin(theta),...
                cos(theta)*sin(phi),...
                cos(phi)*cos(theta),...
                0,...
                0,...
                0,...
                0,...
                0,...
                0,...
                v*cos(phi)*cos(theta) - w*cos(theta)*sin(phi),...
                - u*cos(theta) - w*cos(phi)*sin(theta) - v*sin(phi)*sin(theta),...
                0];

            f10 = [0,...
                0,...
                0,...
                1,...
                sin(phi)*tan(theta),...
                cos(phi)*tan(theta),...
                0,...
                0,...
                0,...
                q*cos(phi)*tan(theta) - r*sin(phi)*tan(theta),...
                r*cos(phi)*(tan(theta)^2 + 1) + q*sin(phi)*(tan(theta)^2 + 1),...
                0];

            f11 = [0,...
                0,...
                0,...
                0,...
                cos(phi),...
                -sin(phi),...
                0,...
                0,...
                0,...
                - r*cos(phi) - q*sin(phi),...
                0,...
                0];

            f12 = [0,...
                0,...
                0,...
                0,...
                sin(phi)/cos(theta),...
                cos(phi)/cos(theta),...
                0,...
                0,...
                0,...
                (q*cos(phi))/cos(theta) - (r*sin(phi))/cos(theta),...
                (r*cos(phi)*sin(theta))/cos(theta)^2 + (q*sin(phi)*sin(theta))/cos(theta)^2,...
                0];
            

            f = [f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11;f12];
            f(13:14,13:14) = [0 1;0 0];

            PHI = expm(f*timeStep);

            cov = PHI*p_m*PHI' + obj.Q;

        end
    end
end

