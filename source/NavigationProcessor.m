classdef NavigationProcessor < handle

    properties
        fullStates;
        engineParameters;
        aircraft;
    end

    properties

        y;
        H;
        K;
        S;
        covariance_plus;
        covariance_minus;
        delStates;
        states_plus;
        R;
        Q;
        

    end

    methods
        function [obj,estPsr,estCarrFreq] = NavigationProcessor(VT,psrRes,carrRes,variances,type)

            % Propagate Satellite Positions
            svStates = VT.svStates;
            if VT.timeSeconds == 0
                obj.fullStates = VT.initialStates;
            end


            switch type
                case 1
                    % Time Update
                    [obj.fullStates,obj.engineParameters] = obj.FVDM_Noise(VT.timeStep,obj.fullStates,obj.engineParameters,VT.controls,obj.aircraft);

                    % Numerically Calculate Jacobian
                    obj.predictedCovariance = obj.calcPredictedCovariance(VT.timeStep);

                    % Predict Covariance

                    % Calculate Estimated (Predicited Pseudorange) and Estimated Carrier Frequency (Predicted Carrier Frequency)



                case 2
                    % Measurement Update
                    obj.y = obj.H*obj.delStates;

                    % Kalman Gain Update
                    obj.S = obj.H*obj.covariance_minus*obj.H' + obj.R;
                    obj.K = obj.covariance_minus*obj.H'*(obj.S)^-1;

                    % Update Estimated States
                    obj.states_plus = obj.K*obj.y;

                    % Covariance Update
                    obj.covariance_plus = (eye(14) - obj.K*obj.H)*obj.covariance_minus;

                    stop = 1;

            end


        end

        function [state,engineParameters] = FVDM_Noise(obj,deltaTime,state,engineParameters,controls,aircraft)
            propR = engineParameters.propR;
            oldFuelFlow = engineParameters.oldFuelFlow;
            oldShaftPower = engineParameters.oldShaftPower;

            Vehicle = aircraft.Vehicle;
            BSFC_LUT = aircraft.BSFC_LUT;
            STGeometry = aircraft.STGeometry;
            refLLA = aircraft.LLA;
            engineForcesVAR = aircraft.engineForcesVAR;
            engineMomentsVAR = aircraft.engineMomentsVAR;
            aeroForcesVAR = aircraft.aeroForcesVAR;
            aeroMomentsVAR = aircraft.aeroMomentsVAR;
            gravityForcesVAR = aircraft.gravityForcesVAR;
            gravityMomentsVAR = aircraft.gravityMomentsVAR;

            % Calculate Atmospheric Parameters for current time step
            LLA = flat2lla(state(7:9)',refLLA(1:2),0,0);
            atmos = Environment('Noise',LLA(1),LLA(2),LLA(3),50,50);

            % Calculate Controller Commands
            latStick = controls(1);
            longStick = controls(2);
            pedalPosn = controls(3);
            throttle = controls(4);


            % Propagate Engine and Prop Forces/Moments
            [~,engineForces,engineMoments,propR,oldFuelFlow,oldShaftPower] = PropagateEngine(atmos,state,throttle,propR,deltaTime,BSFC_LUT,oldFuelFlow,oldShaftPower);
            engineForces = engineForces + randn(1,3)*engineForcesVAR;
            engineMoments = engineMoments + randn(1,3)*engineMomentsVAR;
            % engineForces = 0;
            % engineMoments = 0;

            % Propagate Aerodynamic Forces/Moments
            [~,aeroForces,aeroMoments] = PropagateAero(atmos.density,state,latStick,longStick,pedalPosn,STGeometry);
            aeroForces = aeroForces + randn(1,3)*aeroForcesVAR;
            aeroMoments = aeroMoments + randn(1,3)*aeroMomentsVAR;

            % Propagate Gravity Forces/Moments
            [~,gravityForces,gravityMoments] = PropagateGravity(state,Vehicle.MassProp.Mass,Vehicle.MassProp.r_cg);
            gravityForces = gravityForces + randn(1,3)*gravityForcesVAR;
            gravityMoments = gravityMoments + randn(1,3)*gravityMomentsVAR;

            % Summing All Forces and Moments (Body Frame)
            f_ib_b = engineForces + aeroForces + gravityForces;
            m_ib_b = engineMoments + aeroMoments + gravityMoments;

            % Propgating States
            [~,state,covariance] = PropagateStates(f_ib_b,m_ib_b,Vehicle.MassProp.MOI,Vehicle.MassProp.r_cg,Vehicle.MassProp.InvMassMat,Vehicle.MassProp.Mass,state,deltaTime);

            engineParameters.propR = propR;
            engineParameters.oldFuelFlow = oldFuelFlow;
            engineParameters.oldShaftPower = oldShaftPower;
        end


        function predictedCovariance = calcPredictedCovariance(obj,deltaT,states)

            u = states(1);
            v = states(2);
            w = states(3);
            p = states(4);
            q = states(5);
            r = states(6);
            phi = states(10);
            theta = states(11);
            psi = states(12);

            J(1,1) = (0.023713*p) + (0.0041795*r);
            J(1,2) = 0;
            J(1,3) = 0;
            J(1,4) = -(0.852*q) - (w);
            J(1,5) = (0.0041795*p) - (0.87571*r) + (v);
            J(1,6) = 0;
            J(1,7) = 0;
            J(1,8) = 0;
            J(1,9) = 0;
            J(1,10) = 0;
            J(1,11) = (r);
            J(1,12) = -q;

            J(2,1) = (0.2856*q) + (w);
            J(2,2) = 0;
            J(2,3) = 0;
            J(2,4) = (0.2856*p) - (0.020854*r);
            J(2,5) = -(0.020854*q) - (u);
            J(2,6) = 0;
            J(2,7) = 0;
            J(2,8) = 0;
            J(2,9) = 0;
            J(2,10) = -(r);
            J(2,11) = 0;
            J(2,12) = (p);

            J(3,1) = (0.45051*p) + (0.8609*r) - (v);
            J(3,2) = 0;
            J(3,3) = 0;
            J(3,4) = (0.4*q) + (u);
            J(3,5) = (0.8609*p) - (0.050508*r);
            J(3,6) = 0;
            J(3,7) = 0;
            J(3,8) = 0;
            J(3,9) = 0;
            J(3,10) = (q);
            J(3,11) = -(p);
            J(3,12) = 0;

            J(4,1) = 0;
            J(4,2) = 0;
            J(4,3) = 0;
            J(4,4) = (0.016985*p) + (0.89977*r);
            J(4,5) = (0.89977*q);
            J(4,6) = 0;
            J(4,7) = 0;
            J(4,8) = 0;
            J(4,9) = 0;
            J(4,10) = 0;
            J(4,11) = 0;
            J(4,12) = 0;

            J(5,1) = -(0.11856*p) - (1.0209*r);
            J(5,2) = 0;
            J(5,3) = 0;
            J(5,4) = 0;
            J(5,5) = (0.11856*r) - (1.0209*p);
            J(5,6) = 0;
            J(5,7) = 0;
            J(5,8) = 0;
            J(5,9) = 0;
            J(5,10) = 0;
            J(5,11) = 0;
            J(5,12) = 0;

            J(6,1) = -(0.33754*q);
            J(6,2) = 0;
            J(6,3) = 0;
            J(6,4) = -(0.33754*p) - (0.0018947*r);
            J(6,5) = -(0.0018947*q);
            J(6,6) = 0;
            J(6,7) = 0;
            J(6,8) = 0;
            J(6,9) = 0;
            J(6,10) = 0;
            J(6,11) = 0;
            J(6,12) = 0;

            J(7,1) = 0;
            J(7,2) = w*cos(phi)*sin(psi) + v*sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta) - cos(psi)*sin(phi)*sin(theta);
            J(7,3) =  w*cos(psi)*sin(phi) - v*cos(phi)*cos(psi) - sin(phi)*sin(psi)*sin(theta) - u*cos(theta)*sin(psi) - cos(phi)*sin(psi)*sin(theta);
            J(7,4) = 0;
            J(7,5) = 0;
            J(7,6) = cos(phi)*cos(psi)*cos(theta) - u*cos(psi)*sin(theta) + cos(psi)*cos(theta)*sin(phi);
            J(7,7) = 0;
            J(7,8) = 0;
            J(7,9) = 0;
            J(7,10) = cos(psi)*cos(theta);
            J(7,11) = -cos(phi)*sin(psi);
            J(7,12) = sin(phi)*sin(psi);

            J(8,1) = 0;
            J(8,2) = cos(phi)*sin(ps)*sin(theta) - w*cos(phi)*cos(ps) - v*cos(ps)*sin(phi) - sin(phi)*sin(ps)*sin(theta);
            J(8,3) = u*cos(psi)*cos(theta) - v*cos(phi)*sin(psi) + w*sin(phi)*sin(ps) + cos(phi)*cos(psi)*sin(theta) + cos(psi)*sin(phi)*sin(theta);
            J(8,4) = 0;
            J(8,5) = 0;
            J(8,6) = cos(phi)*cos(theta)*sin(psi) - u*sin(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi);
            J(8,7) = 0;
            J(8,8) = 0;
            J(8,9) = 0;
            J(8,10) = cos(theta)*sin(psi);
            J(8,11) = cos(phi)*cos(psi);
            J(8,12) = -cos(psi)*sin(phi);

            J(9,1) = 0;
            J(9,2) = v*cos(phi)*cos(theta) - w*cos(theta)*sin(phi);
            J(9,3) = 0;
            J(9,4) = 0;
            J(9,5) = 0;
            J(9,6) = - u*cos(theta) - w*cos(phi)*sin(theta) - v*sin(phi)*sin(theta);
            J(9,7) = 0;
            J(9,8) = 0;
            J(9,9) = 0;
            J(9,10) = -sin(theta);
            J(9,11) = cos(theta)*sin(phi);
            J(9,12) = cos(phi)*cos(theta); 

            J(10,1) = 1;
            J(10,2) = q*cos(phi)*tan(theta) - r*sin(phi)*tan(theta);
            J(10,3) = 0;
            J(10,4) = sin(phi)*tan(theta);
            J(10,5) = cos(phi)*tan(theta);
            J(10,6) = r*cos(phi)*(tan(theta)^2 + 1) + q*sin(phi)*(tan(theta)^2 + 1);
            J(10,7) = 0;
            J(10,8) = 0;
            J(10,9) = 0;
            J(10,10) = 0;
            J(10,11) = 0;
            J(10,12) = 0;

            J(11,1) = 0;
            J(11,2) = -r*cos(phi) - q*sin(phi);
            J(11,3) = 0;
            J(11,4) = cos(phi);
            J(11,5) = -sin(phi);
            J(11,6) = 0;
            J(11,7) = 0;
            J(11,8) = 0;
            J(11,9) = 0;
            J(11,10) = 0;
            J(11,11) = 0;
            J(11,12) = 0;

            J(12,1) = 0;
            J(12,2) = (q*cos(phi))/cos(theta) - (r*sin(phi))/cos(theta);
            J(12,3) = 0;
            J(12,4) = sin(phi)/cos(theta);
            J(12,5) = cos(phi)/cos(theta);
            J(12,6) = (r*cos(phi)*sin(theta))/cos(theta)^2 + (q*sin(phi)*sin(theta))/cos(theta)^2;
            J(12,7) = 0;
            J(12,8) = 0;
            J(12,9) = 0;
            J(12,10) = 0;
            J(12,11) = 0;
            J(12,12) = 0;  

            stop = 1;


        end
    end
end

