function [varPsr,varCarr] = calcResVariances(psrRes,R)
pdiTime = 1/50;
chipWidth = 299792458/1.023e6;
cn0 = 10^(45/10);
wavelength = 299792458/1575.42e6;

varPsr = (chipWidth^2/(2*(pdiTime*cn0)^2) + (chipWidth^2.*(psrRes./chipWidth).^2 + 0.25)/(pdiTime*cn0));
varCarr = (wavelength^2/pi/pdiTime)*(2/(pdiTime*cn0)^2 + 2.*R.^2/(pdiTime*cn0));

end