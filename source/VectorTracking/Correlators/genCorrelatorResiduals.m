function [resPsr,resCarr,varPsr,varCarr,newCN0,newAmplitude,newNoise,newPhase] = genCorrelatorResiduals(refPsr,estPsr,refCarr,estCarr,oldCN0,oldAmplitude,oldNoise,activeSVIdx,oldPhase,initCN0)
chipOffset = 0.5;
pdiTime = 1/50;

% Calculate Carrier Frequency and Code Phase Errors
[carrFreqError,codePhaseError] = calcErrors(refPsr,refCarr,estPsr,estCarr);

% Calculate Carrier Phase Errors
phase1 = oldPhase(activeSVIdx);
phase2 = oldPhase(activeSVIdx) + carrFreqError*(pdiTime/2);
phase = (phase1 + phase2)/2;
newPhase = phase2;

% Estimate Amplitude
amplitude = sinc(pi*carrFreqError'*pdiTime);
midAmplitude = sinc(pi*carrFreqError'*pdiTime/2);

% Calculate Noise Estimate
noiseEstimate = sqrt(1./(2*pdiTime*initCN0(activeSVIdx)));

% Calculate Correlators
for i = 1:length(noiseEstimate)
    midIP_codeError(i) = calcCodeError(pdiTime/2,'I',codePhaseError(i),carrFreqError(i),0,phase1(i));
    midQP_codeError(i) = calcCodeError(pdiTime/2,'Q',codePhaseError(i),carrFreqError(i),0,phase1(i));
    IP_codeError(i) = calcCodeError(pdiTime/2,'I',codePhaseError(i),carrFreqError(i),0,phase2(i));
    QP_codeError(i) = calcCodeError(pdiTime/2,'Q',codePhaseError(i),carrFreqError(i),0,phase2(i));
    IE_codeError(i) = calcCodeError(pdiTime,'I',codePhaseError(i),carrFreqError(i),chipOffset,phase(i));
    QE_codeError(i) = calcCodeError(pdiTime,'Q',codePhaseError(i),carrFreqError(i),chipOffset,phase(i));
    IL_codeError(i) = calcCodeError(pdiTime,'I',codePhaseError(i),carrFreqError(i),-chipOffset,phase(i));
    QL_codeError(i) = calcCodeError(pdiTime,'Q',codePhaseError(i),carrFreqError(i),-chipOffset,phase(i));
end

% If Code Error is too big, set to zero
IE_codeError(abs(codePhaseError + chipOffset) > 1) = 0;
QE_codeError(abs(codePhaseError + chipOffset) > 1) = 0;
IL_codeError(abs(codePhaseError - chipOffset) > 1) = 0;
QL_codeError(abs(codePhaseError - chipOffset) > 1) = 0;
midIP_codeError(abs(codePhaseError) > 1) = 0;
midQP_codeError(abs(codePhaseError) > 1) = 0;
IP_codeError(abs(codePhaseError) > 1) = 0;
QP_codeError(abs(codePhaseError) > 1) = 0;

% Generate Correlators
IE = (amplitude.*IE_codeError + noiseEstimate'.*randn(1,length(amplitude)));
QE = (amplitude.*QE_codeError + noiseEstimate'.*randn(1,length(amplitude)));
IL = (amplitude.*IL_codeError + noiseEstimate'.*randn(1,length(amplitude)));
QL = (amplitude.*QL_codeError + noiseEstimate'.*randn(1,length(amplitude)));
midIP = (midAmplitude.*midIP_codeError + sqrt(2).*noiseEstimate'.*randn(1,length(amplitude)));
midQP = (midAmplitude.*midQP_codeError + sqrt(2).* noiseEstimate'.*randn(1,length(amplitude)));
IP = (midAmplitude.*IP_codeError + sqrt(2).*noiseEstimate'.*randn(1,length(amplitude)));
QP = (midAmplitude.*QP_codeError + sqrt(2).*noiseEstimate'.*randn(1,length(amplitude)));


% Generate Discriminators
discFLL = (midIP.*QP - midQP.*IP)/(pi*pdiTime);
discDLL = (IE.^2 + QE.^2 - IL.^2 - QL.^2)/2;

% Estimate Amplitude
power = (IE + IL).^2 + (QE + QL).^2;

% Moving Average Amplitude and Noise Calculation
for i = 1:length(amplitude)
    gpsNoise(i) = var(noiseEstimate(i)*randn(6,1));
end

newAmplitude = 0.99*oldAmplitude(activeSVIdx) + 0.01*power';
newNoise = 0.99*oldNoise(activeSVIdx) + 0.01*gpsNoise';

newCN0 = ((newAmplitude - 4*newNoise)./(2*pdiTime*newNoise));

% Generating pseudorange and pseudorange rate residuals
chipWidth = 299792458/1.023e6;
wavelength = 299792458/1575.42e6;

resPsr = discDLL*chipWidth;
resCarr = discFLL*-wavelength;

% Calculating Variances for Measurement Update
[varPsr,varCarr] = calcResVariances(newCN0);

end

