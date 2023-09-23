function correlator = calcCorrelators(codeError,noiseEstimate,scale,amplitude)

correlator = (amplitude.*codeError + scale*(noiseEstimate)*randn(1));

end