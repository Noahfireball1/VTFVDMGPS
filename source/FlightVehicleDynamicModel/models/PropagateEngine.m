classdef PropagateEngine < handle
    %PROPAGATEENGINE Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = private)

        d = 1.9; % Prop Diameter
        I_prop = 5; % Prop Moment of Inertia

    end

    methods
        function [obj,engineForces,engineMoments,propR,oldFuelFlow,oldShaftPower] = PropagateEngine(atmos,state,throttle,propR,timeStep,BSFC_LUT,oldFuelFlow,oldShaftPower)

            propRPM = propR*(30/pi);

            BSFC = interp1(BSFC_LUT(1,:),BSFC_LUT(2,:),propRPM,"linear");
            if isnan(BSFC)
                BSFC = 0.5695;
            end

            hgPressure = atmos.pressure/3386.39 ;

            airFlowRate = obj.calcAirFlow(atmos.temperature,throttle,hgPressure);
            mixLevel = 0.1;
            [fuelFlow,oldFuelFlow] = obj.calcFuelFlow(mixLevel,airFlowRate,timeStep,oldFuelFlow);

            [delayedShaftPower,oldShaftPower] = obj.calcDelayedShaftPower(fuelFlow,1,BSFC,timeStep,oldShaftPower);

            Q_drive = (delayedShaftPower*746)*propR;
            rho = atmos.density;
            V = norm(state(1:3));
            n = propR/2/pi;

            [engineForces,engineMoments,propRPS] = obj.CTCQPropeff(Q_drive,rho,V,n);

            engineForces = [engineForces;0;0];
            engineMoments = [engineMoments;0;0];
            propR = propRPS*timeStep + propR;

            if propR > 250
                propR = 250;
            elseif propR < 0
                propR = 0;
            end

        end

        function [T,Q,dOmega_dt] = CTCQPropeff(obj,Q_drive,rho,V,n)
            %% Governor State Vector Variable Extraction


            %% Propellor Dynamics (Thrust and Torque)
            J = (V)/(n*obj.d);

            C_T = 0.0041*J^3 - 0.0353*J^2 - 0.0177*J + 0.0947; %thrust coefficient
            C_Q = 0.001*J^4 - 0.008*J^3 + 0.0078*J^2 - 0.0035*J + 0.0114; %torque coefficient

            T = C_T*rho*n^2*obj.d^4;
            Q = C_Q*rho*n^2*obj.d^5;

            %----------RPS-ODE----------
            dOmega_dt = (30/pi)*((Q_drive-Q)/obj.I_prop);

        end

        function airFlowRate = calcAirFlow(~,temp,throttle,hgPressure)

            %% Conversion Factors
            inHg_to_Pa = 3386.39;
            in2_to_m2 = 0.00064516;
            kg_to_lb = 2.20462;
            %% Constants
            gamma = 1.4;
            R = 287;
            AirflowScaler = 0.69 * 0.95;
            DischargeCoeff = 0.6;
            th = 70*throttle;
            D = 2.75;

            Patm_Pa = hgPressure * inHg_to_Pa;

            A_in2 = (pi/4) * D^2 *  (0.01 + (1-0.01)*(1 - cosd(th)));
            A_m2 = A_in2 * in2_to_m2;

            Threshold_X = (2/(gamma+1))^(gamma/(gamma-1));

            x = Threshold_X;

            MAF = sqrt((2*gamma/(gamma-1))*(x^(2/gamma)-x^((gamma+1)/gamma)));

            MassFlux_kgm2s = MAF*Patm_Pa/sqrt(R*temp);
            MassFlow_kgs = MassFlux_kgm2s * A_m2;

            airFlowRate = AirflowScaler*DischargeCoeff*MassFlow_kgs*3600*kg_to_lb;

        end

        function [delayedFuelFlow,oldFuelFlow] = calcFuelFlow(~,mixLevel,airFlowRate,timeStep,oldFuelFlow)

            delayEqn = c2d(tf(1,[0.1 1]),timeStep);

            fuelFlow = mixLevel*airFlowRate;


            delayedFuelFlow = delayEqn.Numerator{1}(2)*fuelFlow*(1/(delayEqn.Denominator{1}(1) - delayEqn.Denominator{1}(2)*oldFuelFlow));
            oldFuelFlow = fuelFlow;
        end

        function [delayedShaftPower,oldShaftPower] = calcDelayedShaftPower(~,fuelFlow,powerFactor,BSFC,timeStep,oldShaftPower)

            delayEqn = c2d(tf(1,[0.05 1]),timeStep);
            shaftPower = powerFactor*fuelFlow/BSFC;

            delayedShaftPower = delayEqn.Numerator{1}(2)*shaftPower*(1/(delayEqn.Denominator{1}(1) - delayEqn.Denominator{1}(2)*oldShaftPower));
            oldShaftPower = shaftPower;
        end


    end
end

