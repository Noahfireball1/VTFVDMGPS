classdef CorrelatorSim < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        cn0 = 10^(45/10);
        pdiTime = 1/50;
        offset = 0.5;
        whiteNoise = randn(1);
        chipWidth = 299792458/1.023e6;
        wavelength = 299792458/1575.42e6;

    end

    methods
        function [obj,psrRes,carrRes,variances] = CorrelatorSim(VT,estPsr,estCarrFreq,sv)

            % Propagate Measurement States
            statesBody = VT.measStates;
            statesECEF = VT.body2ECEF(statesBody);
            statesECEF = [statesECEF 0] + randn(1,7).*[1.5 0.15 1.5 0.15 3.0 0.3 0];

            % Calculate Measured Pseudoranges, Carrier Frequencies
            [measPsr,measCarrFreq] = VT.GE.calcPsr(statesECEF,sv);


            % Calculate Carrier Frequency and Code Phase Errors
            [midCarrFreqError,midCodePhaseError] = obj.calcErrors(measPsr,measCarrFreq,estPsr,estCarrFreq,2);
            [carrFreqError,codePhaseError] = obj.calcErrors(measPsr,measCarrFreq,estPsr,estCarrFreq,1);

            % Calculate Correlators
            whiteNoise = randn(6,1);
            [~,~,midIP,midQP,~,~] = obj.calcCorrelators(midCarrFreqError,midCodePhaseError,whiteNoise);
            [IE,QE,IP,QP,IL,QL] = obj.calcCorrelators(carrFreqError,codePhaseError,whiteNoise);

            % Generate Discriminators
            discFLL = obj.genFLLError(midIP,IP,midQP,QP);
            discPLL = obj.genPLLError(IP,QP);
            discDLL = obj.genDLLError(IE,IL,QE,QL);

            % Estimate Amplitude
            amplitude = (IE + IL).^2 + (QE + QL).^2;

            % Generating pseudorange and pseudorange rate residuals
            psrRes = obj.calcPsrRes(discDLL,amplitude);
            [carrRes,R] = obj.calcCarrRes(discFLL,amplitude,psrRes);

            % Calculating Variances for Measurement Update
            variances = obj.calcResVariances(psrRes,R);



        end

        function [IE,QE,IP,QP,IL,QL] = calcCorrelators(obj,carrFreqError,codePhaseError,whiteNoise)

            noise = sqrt(1/obj.cn0/obj.pdiTime)*whiteNoise;

            IP = sqrt(2*obj.cn0*obj.pdiTime).*sinc((carrFreqError).*obj.pdiTime*2*pi).*(1 - abs(codePhaseError)).*cos(carrFreqError) + noise(1);
            QP = sqrt(2*obj.cn0*obj.pdiTime).*sinc((carrFreqError).*obj.pdiTime*2*pi).*(1 - abs(codePhaseError)).*sin(carrFreqError) + noise(2);
            IE = sqrt(2*obj.cn0*obj.pdiTime).*sinc((carrFreqError).*obj.pdiTime*2*pi).*(1 - abs(codePhaseError - obj.offset)).*cos(carrFreqError) + noise(3);
            QE = sqrt(2*obj.cn0*obj.pdiTime).*sinc((carrFreqError).*obj.pdiTime*2*pi).*(1 - abs(codePhaseError - obj.offset)).*sin(carrFreqError) + noise(4);
            IL = sqrt(2*obj.cn0*obj.pdiTime).*sinc((carrFreqError).*obj.pdiTime*2*pi).*(1 - abs(codePhaseError + obj.offset)).*cos(carrFreqError) + noise(5);
            QL = sqrt(2*obj.cn0*obj.pdiTime).*sinc((carrFreqError).*obj.pdiTime*2*pi).*(1 - abs(codePhaseError + obj.offset)).*sin(carrFreqError) + noise(6);
        end

        function dllError = genDLLError(~,IE, IL, QE, QL)
            dllError = IE.^2 + QE.^2 - IL.^2 - QL.^2;
        end

        function pllError = genPLLError(~,IP, QP)
            pllError = atan(QP./IP)/(2*pi);
        end

        function fllError = genFLLError(obj,lastIP, IP, lastQP, QP)
            cross = lastIP.*QP - IP.*lastQP;
            dot = lastIP.*IP + lastQP.*QP;
            fllError = atan2(cross,dot)./(obj.pdiTime*2*pi);
        end
        function psrRes = calcPsrRes(obj,discDLL,amplitude)
            psrRes = discDLL.*obj.chipWidth./(2.*amplitude);
        end
        function [carrRes,R] = calcCarrRes(obj,discFLL,amplitude,psrRes)
            for i = 1:length(psrRes)
                if psrRes(i)/obj.chipWidth < 1
                    R(i) = 1 - psrRes(i)/obj.chipWidth;
                else
                    R(i) = 0;
                end

                carrRes(i) = discFLL(i)./(-amplitude(i)^2.*R(i)^2*pi*obj.pdiTime*obj.wavelength);
            end
        end
        function variances = calcResVariances(obj,psrRes,R)
            variances.psr = (obj.chipWidth^2/(2*(obj.pdiTime*obj.cn0)^2) + (obj.chipWidth^2.*(psrRes./obj.chipWidth).^2 + 0.25)/(obj.pdiTime*obj.cn0));
            variances.carr = (obj.wavelength^2/pi/obj.pdiTime)*(2/(obj.pdiTime*obj.cn0)^2 + 2.*R.^2/(obj.pdiTime*obj.cn0));

        end
        function [carrFreqError,codePhaseError] = calcErrors(obj,measPsr,measCarrFreq,estPsr,estCarrFreq,type)

            delPsr = measPsr - estPsr;
            delCarr = measCarrFreq - estCarrFreq;

            switch type
                case 1
                    codePhaseError = delPsr/obj.chipWidth;
                    carrFreqError = delCarr/-obj.wavelength;

                case 2
                    codePhaseError = delPsr/obj.chipWidth/2;
                    carrFreqError = delCarr/-obj.wavelength/2;
            end

        end
    end
end