function [carrRes,R] = calcCarrRes(discFLL,amplitude,psrRes)
chipWidth = 299792458/1.023e6;
wavelength = 299792458/1575.42e6;
pdiTime = 1/50;

for i = 1:length(psrRes)
    if psrRes(i)/chipWidth < 1
        R(i) = 1 - psrRes(i)/chipWidth;
    else
        R(i) = 0;
    end

    carrRes(i) = discFLL(i)./(-amplitude(i)^2.*R(i)^2*pi*pdiTime*wavelength);
end
end