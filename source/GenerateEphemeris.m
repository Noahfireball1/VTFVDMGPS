classdef GenerateEphemeris < handle
    %GENERATEEPHEMERIS Summary of this class goes here
    %   Detailed explanation goes here
    properties (Access = public)
        eph;

    end

    properties (Access = private)
        ephYear = 1998;
        ephMonth = 12;
        ephDay = 22;
        login = [];
        netRCFilePath = [];
        c = 299792458;
        wavelength = 299792458/1575.42e6;
        transmitFreq = 1575.42e6;

    end

    properties (Dependent, Access = private)
        julianDay;
        extension;
    end

    methods
        function obj = GenerateEphemeris(date,dir)

            obj.ephYear = year(date);
            obj.ephMonth = month(date);
            obj.ephDay = day(date);

            obj.loadEphemeris(dir)

        end

    end
    methods (Access = private)

        function loadEphemeris(obj,dir)

            % Load or Download Correct Rinex File
            rinexFilePath = obj.loadRinexFile(dir);

            % Parse Ephemeris from File
            switch obj.extension
                case ".rnx.gz"
                    ephTemp = rinexread(rinexFilePath);
                    ephTemp = ephTemp.GPS;

                case ".Z"

            end

            % Select only ephemeris at 0000 hrs
            timeIdxs = ephTemp.Time == datetime(obj.ephYear,obj.ephMonth,obj.ephDay,0,0,0);

            ephTemp = table2struct(ephTemp(timeIdxs,:));
            for i = 1:height(ephTemp)

                obj.eph.satelliteID(1,i) = ephTemp(i).SatelliteID;
                obj.eph.SVClockBias(1,i) = ephTemp(i).SVClockBias;
                obj.eph.SVClockDrift(1,i) = ephTemp(i).SVClockDrift;
                obj.eph.SVClockDriftRate(1,i) = ephTemp(i).SVClockDriftRate;
                obj.eph.IODE(1,i) = ephTemp(i).IODE;
                obj.eph.Crs(1,i) = ephTemp(i).Crs;
                obj.eph.Delta_n(1,i) = ephTemp(i).Delta_n;
                obj.eph.M0(1,i) = ephTemp(i).M0;
                obj.eph.Cuc(1,i) = ephTemp(i).Cuc;
                obj.eph.Eccentricity(1,i) = ephTemp(i).Eccentricity;
                obj.eph.Cus(1,i) = ephTemp(i).Cus;
                obj.eph.sqrtA(1,i) = ephTemp(i).sqrtA;
                obj.eph.Toe(1,i) = ephTemp(i).Toe;
                obj.eph.Cic(1,i) = ephTemp(i).Cic;
                obj.eph.OMEGA0(1,i) = ephTemp(i).OMEGA0;
                obj.eph.Cis(1,i) = ephTemp(i).Cis;
                obj.eph.i0(1,i) = ephTemp(i).i0;
                obj.eph.Crc(1,i) = ephTemp(i).Crc;
                obj.eph.omega(1,i) = ephTemp(i).omega;
                obj.eph.OMEGA_DOT(1,i) = ephTemp(i).OMEGA_DOT;
                obj.eph.IDOT(1,i) = ephTemp(i).IDOT;
                obj.eph.L2ChannelCodes(1,i) = ephTemp(i).L2ChannelCodes;
                obj.eph.GPSWeek(1,i) = ephTemp(i).GPSWeek;
                obj.eph.L2PDataFlag(1,i) = ephTemp(i).L2PDataFlag;
                obj.eph.SVAccuracy(1,i) = ephTemp(i).SVAccuracy;
                obj.eph.SVHealth(1,i) = ephTemp(i).SVHealth;
                obj.eph.TGD(1,i) = ephTemp(i).TGD;
                obj.eph.IODC(1,i) = ephTemp(i).IODC;
                obj.eph.TransmissionTime(1,i) = ephTemp(i).TransmissionTime;
                obj.eph.FitInterval(1,i) = ephTemp(i).FitInterval;
                obj.eph.BRDCOrbit7Spare3(1,i) = ephTemp(i).BRDCOrbit7Spare3;
                obj.eph.BRDCOrbit7Spare4(1,i) = ephTemp(i).BRDCOrbit7Spare4;


            end

        end

        function rinexFile = loadRinexFile(obj,dir)

            yearSuffix = num2str(obj.ephYear);

            baseURL = "https://cddis.nasa.gov/archive/gnss/data/daily/";
            webURL = sprintf('%s/%s/%sn/',num2str(obj.ephYear),string(obj.julianDay),yearSuffix(3:end));
            fileURL = sprintf('AMC400USA_R_%s%s0000_01D_GN.rnx.gz',num2str(obj.ephYear),string(obj.julianDay));
            url = append(baseURL, webURL,fileURL);
            outputFilePath = fullfile(dir.dataGPS,fileURL);
            [path,unzippedFile,~] = fileparts(outputFilePath);

            if ~exist(append(path,filesep,unzippedFile),'file')

                utilities.printText.options(1)

                obj.downloadRinex(url,dir,fileURL)

                obj.uncompressRinex(dir,fileURL)

                utilities.printText.options(2)

            end

            rinexFile = append(path,filesep,unzippedFile);

        end

        function downloadRinex(obj,url,dir,fileURL)

            obj.checkPassword(dir);

            system(sprintf('curl -s -c %scookies.txt --ciphers DEFAULT@SECLEVEL=1 --netrc-file %s -L -O %s',dir.dataGPS,obj.netRCFilePath,url));

            system(sprintf('mv %s %s',fileURL,dir.dataGPS));

            system(sprintf('del %scookies.txt',dir.dataGPS));

        end

        function uncompressRinex(obj,dir,fileURL)

            switch obj.extension
                case ".rnx.gz"
                    gunzip(append(dir.dataGPS,fileURL));
                    system(sprintf('del %s%s',dir.dataGPS,fileURL));
                case ".gz"
                    uncompress(append(dir.dataGPS,fileURL));
            end

        end

        function checkPassword(obj,dir)

            if ~exist(append(dir.config,'CDDISLogin.txt'),"file")

                utilities.printText.options(3)
                username = input('Username:','s');
                password = input('Password:','s');
                obj.createCDDISLoginFile(username,password,dir);
                obj.createNetRCFile(username,password,dir);
            elseif ~exist(append(dir.config,'.netrc'),"file")

                fid = fopen(append(dir.config,'CDDISLogin.txt'));

                loginFile = fscanf(fid,'%c');
                splitLogin = split(loginFile,"password:");

                username = splitLogin{1}(10:end-2);
                password = splitLogin{2}(1:end-2);

                obj.createNetRCFile(username,password,dir);
            end
        end

        function createCDDISLoginFile(~,username,password,dir)

            fid = fopen(append(dir.config,'CDDISLogin.txt'),"wt");
            fprintf(fid,'username:%s\npassword:%s\n',username,password);
            fclose(fid);

        end

        function createNetRCFile(obj,username,password,dir)

            fid = fopen(append(dir.config,'.netrc'),"wt");
            fprintf(fid,'machine urs.earthdata.nasa.gov login %s password %s',username,password);
            fclose(fid);

            obj.netRCFilePath = append(dir.config,'.netrc');

        end

    end

    methods

        function julianDay = get.julianDay(obj)
            julianDay = day(datetime(obj.ephYear,obj.ephMonth,obj.ephDay),"dayofyear");
        end

        function extension = get.extension(obj)
            if obj.ephYear > 2020
                extension = ".rnx.gz";
            else
                extension = ".gz";
            end
        end

    end
    methods (Access = public)
        function [psr,carrFreq] = calcPsr(obj,usrStates, svStates)

            dx = (svStates(1,:) - usrStates(1)).^2;
            dy = (svStates(3,:) - usrStates(3)).^2;
            dz = (svStates(5,:) - usrStates(5)).^2;
            usrVel = [usrStates(2);usrStates(4);usrStates(6)];
            svVel = [svStates(2,:); svStates(4,:); svStates(6,:)];

            range = sqrt(dx + dy + dz);

            psr = range + obj.c*svStates(7,:);
            unitVectors = [sqrt(dx)./range; sqrt(dy)./range; sqrt(dz)./range];

            for i = 1:length(psr)

                carrFreq(i) = obj.transmitFreq - obj.wavelength*(svVel(:,i)'*unitVectors(:,i)) + obj.wavelength*(usrVel'*unitVectors(:,i));

            end

        end

        function [satX,satY,satZ,satU,satV,satW,satClkCorr] = calcSatellitePositions(obj,transmitTime,satelliteID,toe,a_f2,a_f1,a_f0,T_GD,sqrtA,Eccentricity,C_us,C_uc,C_rs,C_rc,C_is,C_ic,M_0,iDOT,i_0,OMEGA,OMEGA_DOT,delta_n,OMEGA_0)

            %% Initialize constants ===================================================
            numOfSatellites = length(satelliteID);

            % GPS constatns

            gpsPi          = 3.1415926535898;  % Pi used in the GPS coordinate
            % system

            %--- Constants for satellite position calculation -------------------------
            Omegae_dot     = 7.2921151467e-5;  % Earth rotation rate, [rad/s]
            GM             = 3.986005e14;      % Earth's universal
            % gravitational parameter,
            % [m^3/s^2]
            F              = -4.442807633e-10; % Constant, [sec/(meter)^(1/2)]

            %% Initialize results =====================================================
            satClkCorr   = zeros(1,numOfSatellites);
            satX = zeros(1,numOfSatellites);
            satY = zeros(1,numOfSatellites);
            satZ = zeros(1,numOfSatellites);
            satU = zeros(1,numOfSatellites);
            satV = zeros(1,numOfSatellites);
            satW = zeros(1,numOfSatellites);

            %% Process each satellite =================================================

            for satNr = 1 : numOfSatellites

                %% Find initial satellite clock correction --------------------------------

                %--- Find time difference ---------------------------------------------
                dt = obj.check_t(transmitTime(satNr) - toe(satNr));

                %--- Calculate clock correction ---------------------------------------
                satClkCorr(satNr) = (a_f2(satNr) * dt + a_f1(satNr)) * dt + ...
                    a_f0(satNr) - ...
                    T_GD(satNr);
                time = transmitTime(satNr) - satClkCorr(satNr);

                %% Find satellite's position ----------------------------------------------

                % Restore semi-major axis
                a   = sqrtA(satNr)*sqrtA(satNr);

                % Time correction
                tk  = obj.check_t(time - toe(satNr));

                % Initial mean motion
                n0  = sqrt(GM/a^3);
                % Mean motion
                n   = n0 + delta_n(satNr);

                % Mean anomaly
                M   = M_0(satNr) + n*tk;
                % Reduce mean anomaly to between 0 and 360 deg
                M   = rem(M + 2*gpsPi, 2*gpsPi);

                % Initial guess of eccentric anomaly
                E   = M;

                %--- Iteratively compute eccentric anomaly ----------------------------
                for ii = 1:10
                    E_old   = E;
                    E       = M + Eccentricity(satNr) * sin(E);
                    dE      = rem(E - E_old, 2*gpsPi);

                    if abs(dE) < 1.e-12
                        % Necessary precision is reached, exit from the loop
                        break;
                    end
                end

                % Reduce eccentric anomaly to between 0 and 360 deg
                E   = rem(E + 2*gpsPi, 2*gpsPi);
                dE = n/(1-Eccentricity(satNr) * cos(E));
                % Compute relativistic correction term
                dtr = F * Eccentricity(satNr) * sqrtA(satNr) * sin(E);

                % Calculate the true anomaly
                nu   = atan2(sqrt(1 - Eccentricity(satNr)^2) * sin(E), cos(E)-Eccentricity(satNr));

                %Compute angle phi
                phi = nu + OMEGA(satNr);
                dphi=sqrt(1-Eccentricity(satNr)^2)*dE/(1-Eccentricity(satNr)*cos(E));
                % Reduce phi to between 0 and 360 deg
                phi = rem(phi, 2*gpsPi);

                % Correct argument of latitude
                u = phi + ...
                    C_uc(satNr) * cos(2*phi) + ...
                    C_us(satNr) * sin(2*phi);
                du=(1+2*(C_us(satNr) * cos(2*phi)-C_uc(satNr) * sin(2*phi)))*dphi;
                % Correct radius
                r = a * (1 - Eccentricity(satNr)*cos(E)) + ...
                    C_rc(satNr) * cos(2*phi) + ...
                    C_rs(satNr) * sin(2*phi);
                dr= a * Eccentricity(satNr) *sin(E) *dE + ...
                    2 * (C_rs(satNr) * cos(2*phi) - C_rc(satNr) * sin(2*phi)) ...
                    * dphi;
                % Correct inclination
                i = i_0(satNr) + iDOT(satNr) * tk + ...
                    C_ic(satNr) * cos(2*phi) + ...
                    C_is(satNr) * sin(2*phi);
                di = 2 * (C_is(satNr) * cos(2*phi) - C_ic(satNr) * sin(2*phi)) ...
                    * dphi + iDOT(satNr);

                % Compute the angle between the ascending node and the Greenwich meridian
                Omega = OMEGA_0(satNr) + (OMEGA_DOT(satNr) - Omegae_dot)*tk - ...
                    Omegae_dot * toe(satNr);
                dOmega = OMEGA_DOT(satNr) - Omegae_dot;
                % Reduce to between 0 and 360 deg
                Omega = rem(Omega + 2*gpsPi, 2*gpsPi);

                %--- Compute satellite coordinates ------------------------------------
                x = cos(u)*r * cos(Omega) - sin(u)*r * cos(i)*sin(Omega);
                y = cos(u)*r * sin(Omega) + sin(u)*r * cos(i)*cos(Omega);
                z = sin(u)*r * sin(i);

                xdash = r * cos(u);
                ydash = r * sin(u);

                dxdash = dr * cos(u) - r * sin(u) * du;
                dydash = dr * sin(u) + r * cos(u) * du;

                Vx = dxdash * cos(Omega) -  dydash * cos(i) * sin(Omega) ...
                    + ydash * sin(Omega) * sin(i) * di - (xdash * sin(Omega) + ...
                    ydash * cos(i) * cos(Omega)) * dOmega;

                Vy = dxdash * sin(Omega) + dydash * cos(i) * cos(Omega) - ...
                    ydash * sin(i) *cos(Omega) * di + (xdash * cos(Omega) - ...
                    ydash * cos(i) * sin(Omega)) * dOmega;

                Vz = dydash * sin(i) + ydash * cos(i) * di;

                satX(satNr) = x;
                satY(satNr) = y;
                satZ(satNr) = z;
                satU(satNr) = Vx;
                satV(satNr) = Vy;
                satW(satNr) = Vz;

                % Include relativistic correction in clock correction --------------------
                satClkCorr(satNr) = (a_f2(satNr) * dt + a_f1(satNr)) * dt + ...
                    a_f0(satNr) - ...
                    T_GD(satNr) + dtr;

            end % for satNr = 1 : numOfSatellites

        end

        function corrTime = check_t(obj,time)

            %CHECKTIME is a subfunction of SATPOS used to check for
            %overlapping time in GPS seconds

            half_week = 302400;     % seconds

            corrTime = time;

            if time > half_week
                corrTime = time - 2*half_week;
            elseif time < -half_week
                corrTime = time + 2*half_week;
            end
        end

    end
end

