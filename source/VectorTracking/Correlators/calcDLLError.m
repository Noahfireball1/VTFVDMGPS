function dllError = calcDLLError(IE, IL, QE, QL)
dllError = IE.^2 + QE.^2 - IL.^2 - QL.^2;
end