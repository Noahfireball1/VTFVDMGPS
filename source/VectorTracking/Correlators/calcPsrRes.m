function psrRes = calcPsrRes(discDLL)
chipWidth = 299792458/1.023e6;

psrRes = discDLL*chipWidth;
end