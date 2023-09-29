function w = formW(S,timestep)

c = 299792458;
h_0 = 2e-25;
h_1 = 7e-25;
h_2 = 6e-25;

clockVar(1,1) = ((h_0/2)*timestep + 2*h_1*timestep^2 + (2/3)*pi^2*h_2*timestep^3);
clockVar(1,2) = (h_1*timestep + pi^2*h_2*timestep^2);
clockVar(2,1) = clockVar(1,2);
clockVar(2,2) = ((h_0/(2*timestep)) + 4*h_1 + (8/3)*pi^2 * h_2*timestep);


w = randn(size(diag(S))).*diag(S);
w(7:8) = sqrt(clockVar)*c*randn(2,1);

end

