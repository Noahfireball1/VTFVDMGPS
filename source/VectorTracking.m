classdef VectorTracking < handle

    properties
        date;
        timeStep;
        startTime;
        endTime;
        timeSeconds;
        svStates;
        rcvrStates;

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

        end

        function [refPsr,refCarrFreq] = process(obj)

            time = obj.startTime:obj.timeStep:obj.endTime;

            for timeIdx = 1:length(time)
                % Prediction and Propagation
                sv = obj.svStates(:,:,timeIdx);
                refStates = obj.rcvrStates(:,timeIdx) + randn(8,1).*[1.5 0.15 1.5 0.15 3.0 0.3 0 0]';

                [refPsr(:,timeIdx),refCarrFreq(:,timeIdx),~] = obj.calcPsr(refStates,sv);

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

        function [psr,carrFreq,unitVectors] = calcPsr(obj,usrStates, svStates)

            dx = (svStates(1,:) - usrStates(1));
            dy = (svStates(3,:) - usrStates(3));
            dz = (svStates(5,:) - usrStates(5));
            usrVel = [usrStates(2);usrStates(4);usrStates(6)];
            svVel = [svStates(2,:); svStates(4,:); svStates(6,:)];

            range = sqrt(dx.^2 + dy.^2 + dz.^2);

            psr = range + obj.c*svStates(7,:);
            unitVectors = [dx./range; dy./range; dz./range];

            for i = 1:length(psr)

                carrFreq(i) = obj.transmitFreq - obj.wavelength*(svVel(:,i)'*unitVectors(:,i)) + obj.wavelength*(usrVel'*unitVectors(:,i));

            end

            % Discard any satellites with a negative elevation
            LLA = ecef2lla([usrStates(1) usrStates(3) usrStates(5)]);
            [~,el,~] = ecef2aer(svStates(1,:),svStates(3,:),svStates(5,:),LLA(1),LLA(2),LLA(3),wgs84Ellipsoid("meter"));

            psr = psr(1,el > 0);
            carrFreq = carrFreq(1,el > 0);
            unitVectors = unitVectors(:,el > 0);

        end

    end

end

