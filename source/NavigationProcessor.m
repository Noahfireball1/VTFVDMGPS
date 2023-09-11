classdef NavigationProcessor < handle

    properties
        fullStates;
    end

    properties

        y;
        H;
        K;
        S;
        covariance_plus;
        covariance_minus = diag([100*ones(1,6),0.00001,0.0000001]);
        delStates;
        states_plus;
        R;
        Q = diag([100*ones(1,12),0.00001,0.0000001]);

    end

    methods
        function [obj,x_m,p_m] = NavigationProcessor(psrRes,carrRes,variances,x_m,p_m,unitVectors)

            % Measurement Update

            infIdx = ~isinf(carrRes);
            carrRes = carrRes(infIdx);
            psrRes = psrRes(infIdx);

            idxOdd = 1:2:2*length(psrRes);
            idxEven = 2:2:2*length(psrRes);
            obj.y(idxOdd) = psrRes;
            obj.y(idxEven) = carrRes;

            [obj.H,obj.R] = obj.formMeasurementMatrix(variances,unitVectors);

            % Kalman Gain Update
            obj.S = obj.H*p_m*obj.H' + obj.R;
            obj.K = p_m*obj.H'*(obj.S)^-1;

            % Update Estimated States
            obj.states_plus = x_m' - obj.K*obj.y';

            % Covariance Update
            obj.covariance_plus = (eye(size(p_m)) - obj.K*obj.H)*p_m;

            p_m = obj.covariance_plus;
            x_m = obj.states_plus;


        end

      
        function [H,R] = formMeasurementMatrix(obj,variances,unitVectors)

            count = 1;
            Rtmp = [];
            for i = 1:2:length(obj.y)
                ux = unitVectors(1,count);
                uy = unitVectors(2,count);
                uz = unitVectors(3,count);

                H(i,1:14) = [ux uy uz 0 0 0 0 0 0 0 0 0 -1 0];
                H(i+1,1:14) = [0 0 0 0 0 0 ux uy uz 0 0 0 0 -1];
                Rtmp = [Rtmp variances.carr(count) variances.psr(count)];

                count = count + 1;          
            end

            R = diag(Rtmp);
        end

    end
end

