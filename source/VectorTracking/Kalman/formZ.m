function Z = formZ(psrRes,carrRes,rcvrStates,estStates)

carrIdx = 1:2:2*length(psrRes);
psrIdx = 2:2:2*length(psrRes);

Z(carrIdx,1) = carrRes;
Z(psrIdx,1) = psrRes;
Z(end+1:end+10,1) = [estStates(4:6) - rcvrStates(4:6);estStates(10:12) - rcvrStates(10:12);estStates(3) - rcvrStates(3);estStates(9) - rcvrStates(9); estStates(1:2) - rcvrStates(1:2)];

end

