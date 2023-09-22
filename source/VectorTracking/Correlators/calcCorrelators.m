function correlator = calcCorrelators(codeError,noiseEstimate,scale,amplitude)

correlator = sqrt(2)*(amplitude.*codeError + scale*noiseEstimate*randn(1,length(amplitude)));

end