function [carrFreqError,codePhaseError] = calcErrors(measPsr,measCarrFreq,estPsr,estCarrFreq)
wavelength = 299792458/1575.42e6;
chipWidth = 299792458/1.023e6;

delPsr = measPsr - estPsr;
delCarr = measCarrFreq - estCarrFreq;


codePhaseError = delPsr/chipWidth;
carrFreqError = delCarr/-wavelength;

end