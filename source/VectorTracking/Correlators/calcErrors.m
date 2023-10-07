function [carrFreqError,codePhaseError] = calcErrors(truePsr,truePsrDot,estPsr,estPsrDot)
wavelength = 299792458/1575.42e6;
chipWidth = 299792458/1.023e6;

delPsr = truePsr - estPsr;
delCarr = truePsrDot - estPsrDot;

codePhaseError = delPsr./chipWidth;
carrFreqError = delCarr./-wavelength;

end