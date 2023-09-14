function Z = formZ(psrRes,carrRes)

idxOdd = 1:2:2*length(psrRes);
idxEven = 2:2:2*length(psrRes);
Z(idxOdd) = carrRes;
Z(idxEven) = psrRes;

end

