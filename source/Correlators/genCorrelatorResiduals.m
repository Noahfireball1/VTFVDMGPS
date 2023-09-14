function [resPsr,resCarr,varPsr,varCarr] = genCorrelatorResiduals(refPsr,estPsr,refCarr,estCarr)

% Calculate Carrier Frequency and Code Phase Errors
[midCarrFreqError,midCodePhaseError] = calcErrors(refPsr,refCarr,estPsr,estCarr,2);
[carrFreqError,codePhaseError] = calcErrors(refPsr,refCarr,estPsr,estCarr,1);

% Calculate Correlators
[~,~,midIP,midQP,~,~] = calcCorrelators(midCarrFreqError,midCodePhaseError);
[IE,QE,IP,QP,IL,QL] = calcCorrelators(carrFreqError,codePhaseError);

% Generate Discriminators
discFLL = calcFLLError(midIP,IP,midQP,QP);
discDLL = calcDLLError(IE,IL,QE,QL);

% Estimate Amplitude
amplitude = (IE + IL).^2 + (QE + QL).^2;

% Generating pseudorange and pseudorange rate residuals
resPsr = calcPsrRes(discDLL,amplitude);
[resCarr,R] = calcCarrRes(discFLL,amplitude,resPsr);

% Calculating Variances for Measurement Update
[varPsr,varCarr] = calcResVariances(resPsr,R);


end

