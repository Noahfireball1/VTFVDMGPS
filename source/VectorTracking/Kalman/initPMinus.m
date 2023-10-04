function pMinus = initPMinus(initialState,clkVar)

phi = initialState(10);
theta = initialState(11);
psi = initialState(12);

cphi = cos(phi);
sphi = sin(phi);

ctheta = cos(theta);
stheta = sin(theta);

cpsi = cos(psi);
spsi = sin(psi);


C_n_b = [1 0 0;0 cphi sphi; 0 -sphi cphi]*...
    [ctheta 0 -stheta; 0 1 0; stheta 0 ctheta]*...
    [cpsi spsi 0; -spsi cpsi 0; 0 0 1];
C_b_n = C_n_b';

pMinus = blkdiag(C_b_n,eye(3),diag([1.5 1.5 3]),C_b_n,clkVar);
end

