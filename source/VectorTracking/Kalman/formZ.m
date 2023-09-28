function Z = formZ(psrRes,carrRes)

idxOdd = 1:2:2*length(psrRes);
idxEven = 2:2:2*length(psrRes);
Z(idxOdd,1) = carrRes;
Z(idxEven,1) = psrRes;

end

