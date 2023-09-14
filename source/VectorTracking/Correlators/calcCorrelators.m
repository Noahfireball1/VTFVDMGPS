function [IE,QE,IP,QP,IL,QL] = calcCorrelators(carrFreqError,codePhaseError)
cn0 = 10^(45/10);
pdiTime = 1/50;
offset = 0.5;

whiteNoise = randn(1,6);
noise = sqrt(1/cn0/pdiTime)*whiteNoise;

IP = sqrt(2* cn0* pdiTime).*sinc((carrFreqError).* pdiTime*2*pi).*(1 - abs(codePhaseError)).*cos(carrFreqError) + noise(1);
QP = sqrt(2* cn0* pdiTime).*sinc((carrFreqError).* pdiTime*2*pi).*(1 - abs(codePhaseError)).*sin(carrFreqError) + noise(2);
IE = sqrt(2* cn0* pdiTime).*sinc((carrFreqError).* pdiTime*2*pi).*(1 - abs(codePhaseError -  offset)).*cos(carrFreqError) + noise(3);
QE = sqrt(2* cn0* pdiTime).*sinc((carrFreqError).* pdiTime*2*pi).*(1 - abs(codePhaseError -  offset)).*sin(carrFreqError) + noise(4);
IL = sqrt(2* cn0* pdiTime).*sinc((carrFreqError).* pdiTime*2*pi).*(1 - abs(codePhaseError +  offset)).*cos(carrFreqError) + noise(5);
QL = sqrt(2* cn0* pdiTime).*sinc((carrFreqError).* pdiTime*2*pi).*(1 - abs(codePhaseError +  offset)).*sin(carrFreqError) + noise(6);
end