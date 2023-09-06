classdef SatellitePositions < handle

    methods
        function [obj,svStates] = SatellitePositions(VT,rinex)

            timeArray = 0:VT.timeStep:VT.endTime;

            dayVector = repmat(day(VT.date),length(timeArray),1);
            monthVector = repmat(month(VT.date),length(timeArray),1);
            yearVector = repmat(year(VT.date),length(timeArray),1);
            hourVector = zeros(length(timeArray),1);
            minuteVector = zeros(length(timeArray),1);
            secondVector = timeArray';

            dateVector = [yearVector monthVector dayVector hourVector minuteVector secondVector];

            dateTimeVector = datetime(dateVector,'Format','d-MMM-y HH:mm:ss.SSS');
            data = rinex.rinex.GPS;
            [~,satIdx] = unique(data.SatelliteID);
            data = data(satIdx,:);

            utilities.printText.options(6);
            upd = utilities.progressbar.textprogressbar(length(timeArray),'startmsg','');
            for i = 1:length(timeArray)
                upd(i)
                [satPos,satVel,satID] = gnssconstellation(dateTimeVector(i),data,"GNSSFileType","RINEX");

                svStates(:,:,i) = [satPos(:,1)';satVel(:,1)';satPos(:,2)';satVel(:,2)';satPos(:,3)';satVel(:,3)';zeros(1,length(satID))];
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

