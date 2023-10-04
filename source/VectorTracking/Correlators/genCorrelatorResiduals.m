function [resPsr,resCarr,varPsr,varCarr,newCN0,newAmplitude,newNoise,newPhase] = genCorrelatorResiduals(refPsr,estPsr,refCarr,estCarr,oldCN0,oldAmplitude,oldNoise,activeSVIdx,oldPhase,initCN0)
chipOffset = 0.5;
pdiTime = 1/50;

% Calculate Carrier Frequency and Code Phase Errors
[carrFreqError,codePhaseError] = calcErrors(refPsr,refCarr,estPsr,estCarr);

% Calculate Carrier Phase Errors
phase1 = oldPhase(activeSVIdx)';
phase2 = oldPhase(activeSVIdx)' + carrFreqError*(pdiTime/2);
phase = (phase1 + phase2)/2;
newPhase = phase2;

% Estimate Amplitude
amplitude = sinc(pi*carrFreqError*pdiTime);
midAmplitude = sinc(pi*carrFreqError*pdiTime/2);


% Calculate Noise Estimate
noiseEstimate = sqrt(1./(2*pdiTime*initCN0(activeSVIdx)));

% Calculate Correlators
for i = 1:length(noiseEstimate)
    midIP(i) = calcCodeError(pdiTime/2,'I',codePhaseError(i),carrFreqError(i),0,phase2(i),noiseEstimate(i),midAmplitude(i),sqrt(2));
    midQP(i) = calcCodeError(pdiTime/2,'Q',codePhaseError(i),carrFreqError(i),0,phase2(i),noiseEstimate(i),midAmplitude(i),sqrt(2));
    IP(i) = calcCodeError(pdiTime/2,'I',codePhaseError(i),carrFreqError(i),0,phase1(i),noiseEstimate(i),midAmplitude(i),sqrt(2));
    QP(i) = calcCodeError(pdiTime/2,'Q',codePhaseError(i),carrFreqError(i),0,phase1(i),noiseEstimate(i),midAmplitude(i),sqrt(2));
    IE(i) = calcCodeError(pdiTime,'I',codePhaseError(i),carrFreqError(i),-chipOffset,phase(i),noiseEstimate(i),amplitude(i),1);
    QE(i) = calcCodeError(pdiTime,'Q',codePhaseError(i),carrFreqError(i),-chipOffset,phase(i),noiseEstimate(i),amplitude(i),1);
    IL(i) = calcCodeError(pdiTime,'I',codePhaseError(i),carrFreqError(i),chipOffset,phase(i),noiseEstimate(i),amplitude(i),1);
    QL(i) = calcCodeError(pdiTime,'Q',codePhaseError(i),carrFreqError(i),chipOffset,phase(i),noiseEstimate(i),amplitude(i),1);
end

% Generate Discriminators
discFLL = calcFLLError(midIP,IP,midQP,QP);
discDLL = calcDLLError(IE,IL,QE,QL);

% Estimate Amplitude
power = ((IE + IL).^2 + (QE + QL).^2);

% Moving Average Amplitude and Noise Calculation
for i = 1:length(amplitude)
    gpsNoise(i) = var(noiseEstimate(i)*randn(6,1));
end

newAmplitude = 0.99*oldAmplitude(activeSVIdx)' + 0.01*power;
newNoise = 0.99*oldNoise(activeSVIdx)' + 0.01*gpsNoise;

newCN0 = ((newAmplitude - 4*newNoise)./(2*pdiTime*newNoise));

% Generating pseudorange and pseudorange rate residuals
resPsr = calcPsrRes(discDLL);
resCarr = calcCarrRes(discFLL);

% Calculating Variances for Measurement Update
[varPsr,varCarr] = calcResVariances(newCN0);

end

