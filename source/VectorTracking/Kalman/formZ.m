function Z = formZ(psrRes,carrRes,rcvrStates)

carrIdx = 1:2:2*length(psrRes);
psrIdx = 2:2:2*length(psrRes);

Z(carrIdx,1) = carrRes;
Z(psrIdx,1) = psrRes;

end

