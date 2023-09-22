function [carrRes] = calcCarrRes(discFLL)

wavelength = 299792458/1575.42e6;

carrRes = discFLL*-wavelength;

end