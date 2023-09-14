function [carrFreqError,codePhaseError] = calcErrors(measPsr,measCarrFreq,estPsr,estCarrFreq,type)
wavelength = 299792458/1575.42e6;
chipWidth = 299792458/1.023e6;

delPsr = measPsr - estPsr;
delCarr = measCarrFreq - estCarrFreq;

switch type
    case 1
        codePhaseError = delPsr/chipWidth;
        carrFreqError = delCarr/-wavelength;

    case 2
        codePhaseError = delPsr/chipWidth/2;
        carrFreqError = delCarr/-wavelength/2;
end

end