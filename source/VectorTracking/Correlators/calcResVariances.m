function [varPsr,varCarr] = calcResVariances(cn0)
pdiTime = 1/50;
chipWidth = 299792458/1.023e6;
wavelength = 299792458/1575.42e6;


varPsr = chipWidth^2*(1./(2*pdiTime^2*cn0.^2) + 1./(4*pdiTime*cn0));

varCarr = (wavelength/(pi*pdiTime))^2 ...
*(2./(pdiTime^2*cn0.^2) + 2./(pdiTime*cn0));
                                                    

end