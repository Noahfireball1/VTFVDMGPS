function Qd = formQ(states,timestep)
c = 299792458;
h_0 = 2e-25;
h_1 = 7e-25;
h_2 = 6e-25;

Qd = zeros(14,14);

Qd(13,13) = ((h_0/2)*timestep + 2*h_1*timestep^2 + (2/3)*pi^2*h_2*timestep^3)*c;
Qd(13,14) = (h_1*timestep + pi^2*h_2*timestep^2)*c;
Qd(14,13) = Qd(13,14);
Qd(14,14) = ((h_0/(2*timestep)) + 4*h_1 + (8/3)*pi^2 * h_2*timestep)*c;
end

