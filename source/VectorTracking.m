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
        unitVectors;
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

            start = datetime(2022,12,22,0,0,0,0);
            stop = datetime(2022,12,22,0,6,30,0);

            sc = satelliteScenario(start,stop,obj.timeStep,"AutoSimulate",false);
            sat = satellite(sc,rinexread(obj.GE.rinexFilePath),"OrbitPropagator","gps");
            [pos,vel] = states(sat,start,"CoordinateFrame","ecef")

            obj.ephemeris = obj.GE.eph;

        end

        function process(obj)

            % Initialize Variables
            time = obj.startTime:obj.timeStep:obj.endTime;
            ecef = lla2ecef([32.6099 85.4808 232],'WGS84');
            X_m = [ecef(1) 0 ecef(2) 0 ecef(3) 0 0 0]';
            P_m = diag([1.5 0.15 1.5 0.15 3.0 0.30 0 0]);
            sv = obj.calcSVStates(0);
            sv = obj.sortSVs(X_m,sv);
            [estPsr,estCarrFreq,obj.unitVectors] = obj.GE.calcPsr(X_m,sv);
            refPsr = estPsr;
            refCarrFreq = estCarrFreq;

            VLL = VDFLL();

            transmitTime = 0;
            for timeIdx = 1:length(time)

                %% Simulating Correlators
                [CS,psrRes,carrRes,variances] = CorrelatorSim(estPsr,estCarrFreq,refPsr,refCarrFreq);

                %% Navigation Processor
                % Time Update
                [X_m,P_m] = VLL.timeUpdate(1/50,X_m,P_m);

                % Measurement Update
                [X_p,P_p] = VLL.measurementUpdate(obj,X_m,P_m,psrRes,carrRes,variances);
                X_m = X_p;
                P_m = P_p;

                % Prediction and Propagation
                transmitTime = transmitTime + obj.timeStep;
                sv = obj.calcSVStates(transmitTime);
                sv = obj.sortSVs(X_m,sv);
                [estPsr,estCarrFreq,obj.unitVectors] = obj.GE.calcPsr(X_m,sv);

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

        function [psr,carrFreq,unitVectors] = initialization(obj)

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

            sv = [satX;satU;satY;satV;satZ;satW;satClkCorr];

            ecefStates = obj.body2ECEF(obj.initialStates);

            sv = obj.sortSVs(ecefStates,sv);

            [psr,carrFreq,unitVectors] = obj.GE.calcPsr(rcvrStates,sv);

        end

        function svStates = calcSVStates(obj,navTime)

            transmitTime = obj.ephemeris.Toe + navTime;
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

        function svs = sortSVs(obj,states,svs)

            svX = svs(1,:);
            svY = svs(3,:);
            svZ = svs(5,:);

            [LLA] = ecef2lla([states(1) states(3) states(5)]);

            [~,el,~] = ecef2aer(svX,svY,svZ,LLA(1),LLA(2),LLA(3),wgs84Ellipsoid("meter"));

            svs = svs(:,el > 10);
        end

    end

end

