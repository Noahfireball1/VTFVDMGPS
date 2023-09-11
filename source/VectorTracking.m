classdef VectorTracking < handle

    properties
        date;
        timeStep;
        startTime;
        endTime;
        timeSeconds;
        initialStates;
        svStates;
        rcvrStates;
        controls;
        RP;

    end

    properties (Constant, Access = private)
        c = 299792458;
        transmitFreq = 1575.42e6;
        wavelength = 299792458/1575.42e6;

    end

    methods
        function obj = VectorTracking(configFilePath,dirs)

            % Read Selected Configuration File
            config = utilities.yamlmatlab.ReadYaml(configFilePath);

            % Parse Yaml File
            obj.parseYaml(config,dirs);

            % Download Ephemeris File from Internet
            rinex = GenerateEphemeris(obj.date,dirs);

            % Generate Satellite Positions for Simulation Time
            [~,obj.svStates] = SatellitePositions(obj,rinex);

            % Generate Receiever Positions for Simulation Time
            load("rcvr_Scurve.mat")
            obj.rcvrStates = rcvrStates;
            obj.controls = controls;

        end

        function process(obj)

            time = obj.startTime:obj.timeStep:obj.endTime;
            measTime = 0:1/50:time(end);
            measIdx = 1;
            x_m_N = obj.initialStates;
            p_m = zeros(14);

            for timeIdx = 1:length(time)

                % Time Update
                cont = obj.controls(:,measIdx);
                bodyStates = obj.RP.noiseFVDM(x_m_N,cont);

                x_m = obj.RP.body2ecef(bodyStates);
                p_m = obj.RP.covariance(obj.timeStep,bodyStates,p_m);

                % Measurement Update
                if (measTime(measIdx) - time(timeIdx)) == 0

                    % Prediction and Propagation
                    sv = obj.svStates(:,:,measIdx);
                    refStates = obj.rcvrStates(:,measIdx) + randn(8,1).*[0.15 0.15 0.30 1.5 1.5 3.0 0 0]';

                    [refPsr,refCarrFreq,uv] = obj.calcPsr(refStates,sv);
                    [estPsr,estCarrFreq] = obj.calcPsr([x_m(1:3) x_m(7:9)],sv);

                    % Simulating Correlators
                    [~,psrRes,carrRes,variances] = CorrelatorSim(estPsr,estCarrFreq,refPsr,refCarrFreq);

                    % Extended Kalman Filter
                    [~,x_m,p_m] = NavigationProcessor(psrRes,carrRes,variances,x_m,p_m,uv);

                    measIdx = measIdx + 1;
                end
                x_m_N = obj.RP.ecef2body(x_m);


            end

        end
    end

    methods (Access = private)
        function parseYaml(obj,config,dirs)

            gen = config.general;

            obj.timeStep = 1/config.aircraft.frequency;
            obj.startTime = 0;
            obj.endTime = config.aircraft.time;
            obj.date = datetime(gen.year,gen.month,gen.day);

            obj.initialStates = [config.aircraft.initialState.u;...
                config.aircraft.initialState.w;...
                config.aircraft.initialState.v;...
                config.aircraft.initialState.p;...
                config.aircraft.initialState.q;...
                config.aircraft.initialState.r;...
                config.aircraft.initialState.x;...
                config.aircraft.initialState.y;...
                config.aircraft.initialState.z;...
                config.aircraft.initialState.phi;...
                str2num(config.aircraft.initialState.theta);...
                str2num(config.aircraft.initialState.psi)];

            Vehicle = load("DA40.mat");
            BSFC_LUT = load("DA40ENGINE.mat");
            STGeometry = load("DA40STGEOM.mat");
            selWaypoints = load(sprintf('%s.mat',config.aircraft.waypoints));

            aircraft.LLA = [selWaypoints.refLL -obj.initialStates(9)];
            aircraft.BSFC_LUT = BSFC_LUT.BSFC_LUT;
            aircraft.Vehicle = Vehicle.Vehicle;
            aircraft.STGeometry = STGeometry.ST_Geometry;
            aircraft.lookaheadDist = config.aircraft.lookaheadDistance;
            aircraft.waypoints = selWaypoints.waypoints;
            aircraft.engineForcesVAR = config.noise.engineForcesVAR;
            aircraft.engineMomentsVAR = config.noise.engineMomentsVAR;
            aircraft.aeroForcesVAR = config.noise.aeroForcesVAR;
            aircraft.aeroMomentsVAR = config.noise.aeroMomentsVAR;
            aircraft.gravityForcesVAR = config.noise.gravityForcesVAR;
            aircraft.gravityMomentsVAR = config.noise.gravityMomentsVAR;

            obj.RP = ReceiverPositions(obj.startTime,obj.endTime,obj.timeStep,aircraft);


        end

    end

    methods (Access = public)


        function [psr,carrFreq,unitVectors] = calcPsr(obj,usrStates, svStates)

            dx = (svStates(2,:) - usrStates(4));
            dy = (svStates(4,:) - usrStates(5));
            dz = (svStates(6,:) - usrStates(6));
            usrVel = [usrStates(1);usrStates(2);usrStates(3)];
            svVel = [svStates(2,:); svStates(4,:); svStates(6,:)];

            range = sqrt(dx.^2 + dy.^2 + dz.^2);

            psr = range + obj.c*svStates(7,:);
            unitVectors = [dx./range; dy./range; dz./range];

            for i = 1:length(psr)

                carrFreq(i) = obj.transmitFreq - obj.wavelength*(svVel(:,i)'*unitVectors(:,i)) + obj.wavelength*(usrVel'*unitVectors(:,i));

            end

            % Discard any satellites with a negative elevation
            LLA = ecef2lla([usrStates(4) usrStates(5) usrStates(6)]);
            [~,el,~] = ecef2aer(svStates(1,:),svStates(3,:),svStates(5,:),LLA(1),LLA(2),LLA(3),wgs84Ellipsoid("meter"));

            psr = psr(1,el > 0);
            carrFreq = carrFreq(1,el > 0);
            unitVectors = unitVectors(:,el > 0);

        end

    end

end

