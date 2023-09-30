function Qd = formQ(S,G,timestep)
c = 299792458;
h_0 = 2e-25;
h_1 = 7e-25;
h_2 = 6e-25;

clkVar(1,1) = ((h_0/2)*timestep + 2*h_1*timestep^2 + (2/3)*pi^2*h_2*timestep^3);
clkVar(1,2) = (h_1*timestep + pi^2*h_2*timestep^2);
clkVar(2,1) = clkVar(1,2);
clkVar(2,2) = ((h_0/(2*timestep)) + 4*h_1 + (8/3)*pi^2 * h_2*timestep);

SwClk = blkdiag(sqrt(S)./timestep,clkVar*c*c*0);

Qd = (G*SwClk*G')*timestep;
end
