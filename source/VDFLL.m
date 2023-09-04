classdef VDFLL < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        A_d;
        P_m = zeros(8);
        Q_d;
        H;
        L;
        P_p;
        X_m;
        X_p;
        R_d;
        y;
        S;
        h_0 = 2e-25;
        h_1 = 7e-25;
        h_2 = 6e-25;
    end

    properties (Dependent)
        S_f;
        S_g;

    end

    properties (Constant)
        c = 299792458;

    end

    methods (Access = public)
        function obj = VDFLL()

        end

        function [X_m,P_m] = timeUpdate(obj,timeStep,X_m,P_m)

            k_d = [1 timeStep;0 1];
            obj.A_d = blkdiag(k_d,k_d,k_d,k_d);
            q_d = [timeStep^3/3 timeStep^2/2; timeStep^2/2 timeStep];
            obj.Q_d = blkdiag(q_d,q_d,q_d,...
                [obj.S_f*timeStep + obj.S_g*timeStep^3/3, obj.S_g*timeStep^2/2;obj.S_g*timeStep^2/2, obj.S_g*timeStep]);
            X_m = obj.A_d*X_m;
            P_m = obj.A_d*P_m*obj.A_d' + obj.Q_d;

        end

        function [X_p,P_p] = measurementUpdate(obj,VT,X_m,P_m,psrRes,carrRes,variances)

            infIdx = ~isinf(carrRes);
            carrRes = carrRes(infIdx);
            psrRes = psrRes(infIdx);

            variances.psr = variances.psr(infIdx);
            variances.carr = variances.carr(infIdx);
            unitVectors = VT.unitVectors(:,infIdx);

            obj.formMeasurementMatrix(variances,unitVectors);

            idxOdd = 1:2:2*length(psrRes);
            idxEven = 2:2:2*length(psrRes);
            obj.y(idxOdd) = psrRes;
            obj.y(idxEven) = carrRes;

            obj.S = obj.H*P_m*obj.H' + obj.R_d;
            obj.L = P_m*obj.H'*(obj.S)^-1;

            % Update Estimated States
            X_p = X_m - obj.L*obj.y';
            P_p = (eye(size(P_m)) - obj.L*obj.H)*P_m;

        end

    end

    methods (Access = private)

        function formMeasurementMatrix(obj,variances,unitVectors)

            count = 1;
            for i = 1:2:2*width(variances.psr)
                ux = unitVectors(1,count);
                uy = unitVectors(2,count);
                uz = unitVectors(3,count);

                obj.H(i,:) = [ux 0 uy 0 uz 0 -1 0];
                obj.H(i+1,:) = [0 ux 0 uy 0 uz 0 -1];
                obj.R_d(i,i) = variances.psr(count);
                obj.R_d(i+1,i+1) = variances.carr(count);

                count = count + 1;
            end

        end

    end

    methods
        function S_f = get.S_f(obj)

            S_f = obj.c^2*obj.h_0/2;

        end

        function S_g = get.S_g(obj)

            S_g = obj.c^2*2*pi^2*obj.h_2;

        end

    end
end