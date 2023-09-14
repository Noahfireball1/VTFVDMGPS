function psrRes = calcPsrRes(discDLL,amplitude)
chipWidth = 299792458/1.023e6;

psrRes = discDLL.*chipWidth./(2.*amplitude);
end