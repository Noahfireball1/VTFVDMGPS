function dllError = calcDLLError(IE, IL, QE, QL)
E = sqrt(IE.^2 + QE.^2);
L = sqrt(IL.^2 + QL.^2);

dllError = 0.5*((E-L)/(E+L));
end