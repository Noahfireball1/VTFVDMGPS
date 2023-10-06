function [carrFreqError,codePhaseError] = calcErrors(truePsr,truePsrDot,estPsr,estPsrDot)
wavelength = 299792458/1575.42e6;
chipWidth = 299792458/1.023e6;

trueCarrier = (truePsrDot/-wavelength + 1575.42e6);
trueChipWidth = 299792458./((1.023e6/1575.42e6)*trueCarrier);
trueWaveLength = 299792458./trueCarrier;


delPsr = truePsr - estPsr;
delCarr = truePsrDot - estPsrDot;

codePhaseError = delPsr./trueChipWidth;
carrFreqError = delCarr./-trueWaveLength;

hold on
plot(trueChipWidth)

end