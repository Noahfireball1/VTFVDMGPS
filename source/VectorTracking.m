classdef VectorTracking < handle

    properties
        date;
        timeStep;
        startTime;
        endTime;
        timeSeconds;
        ephemeris;
        aircraft;
        engineParameters;
        initialStates;
        measStates;
        controls;
        svStates;
        GE;


    end

    properties (Access = private)



    end

    methods
        function obj = VectorTracking(configFilePath,dirs)

            % Read Selected Configuration File
            config = utilities.yamlmatlab.ReadYaml(configFilePath);

            % Parse Yaml File
            obj.parseYaml(config,dirs);

            % Download Ephemeris File from Internet
            obj.GE = GenerateEphemeris(obj.date,dirs);
            obj.ephemeris = obj.GE.eph;

        end

        function process(obj)

            % Initialize Variables
            time = obj.startTime:obj.timeStep:obj.endTime;
            measTime = obj.startTime:1/50:obj.endTime;
            [estPsr,estCarrFreq] = obj.initialization;

            WaypointFollower = uavWaypointFollower('UAVType','fixed-wing','StartFrom','first','Waypoints',flipud(obj.aircraft.waypoints),'TransitionRadius',1000);
            state = obj.initialStates;

            measIdx = 1;
            for timeIdx = 1:length(time)

                obj.timeSeconds = time(timeIdx);
                obj.svStates = obj.calcSVStates();

                if (measTime(measIdx) - time(timeIdx)) == 0
                    %% Generating Measurement States and Controls
                    [state,obj.engineParameters,trueControls] = FVDM_Truth(obj.timeStep,state,obj.engineParameters,WaypointFollower,obj.aircraft);
                    obj.measStates = state;
                    obj.controls = trueControls;
                    
                    %% Simulating Correlators
                    [CS,psrRes,carrRes,variances] = CorrelatorSim(obj,estPsr,estCarrFreq);

                    measIdx = measIdx + 1;

                    [NP,estPsr,estCarrFreq] = NavigationProcessor(obj,psrRes,carrRes,variances,2);
                else
                    %% Navigation Processor
                    [NP,estPsr,estCarrFreq] = NavigationProcessor(obj,psrRes,carrRes,variances,1);
                end
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
                config.aircraft.initialState.v;...
                config.aircraft.initialState.w;...
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

            obj.aircraft.LLA = [selWaypoints.refLL -obj.initialStates(9)];
            obj.aircraft.BSFC_LUT = BSFC_LUT.BSFC_LUT;
            obj.aircraft.Vehicle = Vehicle.Vehicle;
            obj.aircraft.STGeometry = STGeometry.ST_Geometry;
            obj.aircraft.lookaheadDist = config.aircraft.lookaheadDistance;
            obj.aircraft.waypoints = selWaypoints.waypoints;
            obj.aircraft.engineForcesVAR = config.noise.engineForcesVAR;
            obj.aircraft.engineMomentsVAR = config.noise.engineMomentsVAR;
            obj.aircraft.aeroForcesVAR = config.noise.aeroForcesVAR;
            obj.aircraft.aeroMomentsVAR = config.noise.aeroMomentsVAR;
            obj.aircraft.gravityForcesVAR = config.noise.gravityForcesVAR;
            obj.aircraft.gravityMomentsVAR = config.noise.gravityMomentsVAR;

            obj.engineParameters.oldFuelFlow = 0;
            obj.engineParameters.oldShaftPower = 0;
            obj.engineParameters.propR = 200;


        end

        function [psr,carrFreq] = initialization(obj)

            body2ned = genDCM('rad',obj.initialStates(10:12),[3 2 1]);
            nedP = body2ned*obj.initialStates(7:9);
            nedV = body2ned*obj.initialStates(1:3);
            clockBias = 0;
            clockDrift = 0;

            [x,y,z] = ned2ecef(nedP(1),nedP(2),nedP(3),obj.aircraft.LLA(1),obj.aircraft.LLA(2),obj.aircraft.LLA(3),wgs84Ellipsoid("meter"));
            [xdot,ydot,zdot] = ned2ecefv(nedV(1),nedV(2),nedV(3),obj.aircraft.LLA(1),obj.aircraft.LLA(2));

            rcvrStates = [x;xdot;y;ydot;z;zdot;clockBias;clockDrift] + [randn(6,1).*[1.5;1.5;3;0.15;0.15;0.3];0;0];

            transmitTime = obj.ephemeris.Toe;
            satelliteID = obj.ephemeris.satelliteID;
            toe = obj.ephemeris.Toe;
            a_f2 = obj.ephemeris.SVClockDriftRate;
            a_f1 = obj.ephemeris.SVClockDrift;
            a_f0 = obj.ephemeris.SVClockDrift;
            T_GD = obj.ephemeris.TGD;
            sqrtA = obj.ephemeris.sqrtA;
            Eccentricity = obj.ephemeris.Eccentricity;
            C_us = obj.ephemeris.Cus;
            C_uc = obj.ephemeris.Cuc;
            C_rs = obj.ephemeris.Crs;
            C_rc = obj.ephemeris.Crc;
            C_is = obj.ephemeris.Cis;
            C_ic = obj.ephemeris.Cic;
            M_0 = obj.ephemeris.M0;
            iDOT = obj.ephemeris.IDOT;
            i_0 = obj.ephemeris.i0;
            OMEGA = obj.ephemeris.omega;
            OMEGA_DOT = obj.ephemeris.OMEGA_DOT;
            delta_n = obj.ephemeris.Delta_n;
            OMEGA_0 = obj.ephemeris.OMEGA0;

            [satX,satY,satZ,satU,satV,satW,satClkCorr] = obj.GE.calcSatellitePositions(transmitTime,satelliteID,toe,a_f2,a_f1,a_f0,T_GD,sqrtA,Eccentricity,C_us,C_uc,C_rs,C_rc,C_is,C_ic,M_0,iDOT,i_0,OMEGA,OMEGA_DOT,delta_n,OMEGA_0);

            svStates = [satX;satU;satY;satV;satZ;satW;satClkCorr];

            [psr,carrFreq] = obj.GE.calcPsr(rcvrStates,svStates);

        end

        function svStates = calcSVStates(obj)

            transmitTime = obj.ephemeris.Toe + obj.timeSeconds;
            satelliteID = obj.ephemeris.satelliteID;
            toe = obj.ephemeris.Toe;
            a_f2 = obj.ephemeris.SVClockDriftRate;
            a_f1 = obj.ephemeris.SVClockDrift;
            a_f0 = obj.ephemeris.SVClockDrift;
            T_GD = obj.ephemeris.TGD;
            sqrtA = obj.ephemeris.sqrtA;
            Eccentricity = obj.ephemeris.Eccentricity;
            C_us = obj.ephemeris.Cus;
            C_uc = obj.ephemeris.Cuc;
            C_rs = obj.ephemeris.Crs;
            C_rc = obj.ephemeris.Crc;
            C_is = obj.ephemeris.Cis;
            C_ic = obj.ephemeris.Cic;
            M_0 = obj.ephemeris.M0;
            iDOT = obj.ephemeris.IDOT;
            i_0 = obj.ephemeris.i0;
            OMEGA = obj.ephemeris.omega;
            OMEGA_DOT = obj.ephemeris.OMEGA_DOT;
            delta_n = obj.ephemeris.Delta_n;
            OMEGA_0 = obj.ephemeris.OMEGA0;

            [satX,satY,satZ,satU,satV,satW,satClkCorr] = obj.GE.calcSatellitePositions(transmitTime,satelliteID,toe,a_f2,a_f1,a_f0,T_GD,sqrtA,Eccentricity,C_us,C_uc,C_rs,C_rc,C_is,C_ic,M_0,iDOT,i_0,OMEGA,OMEGA_DOT,delta_n,OMEGA_0);

            svStates = [satX;satU;satY;satV;satZ;satW;satClkCorr];

        end
    end

    methods (Access = public)
        function states_ECEF = body2ECEF(obj,state)
            DCM_b_n = genDCM('rad',state(10:12),[3 2 1]);

            pos_NED = DCM_b_n*state(7:9);
            vel_NED = DCM_b_n*state(1:3);

            [x,y,z] = ned2ecef(pos_NED(1),pos_NED(2),pos_NED(3),obj.aircraft.LLA(1),obj.aircraft.LLA(2),0,wgs84Ellipsoid("meter"),"degrees");
            [u,v,w] = ned2ecefv(vel_NED(1),vel_NED(2),vel_NED(3),obj.aircraft.LLA(1),obj.aircraft.LLA(2),"degrees");

            states_ECEF = [x,u,y,v,z,w];
        end

    end

end

