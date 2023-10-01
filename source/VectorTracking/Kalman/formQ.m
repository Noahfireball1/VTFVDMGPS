function Qd = formQ(stateVariance,clockType,timestep)

c = 299792458;
switch upper(clockType)
    case 'LOW-TCXO'
        h_0 = 2e-19;
        h_1 = 7e-21;
        h_2 = 2e-20;

    case 'HIGH-TCXO'
        h_0 = 2e-21;
        h_1 = 1e-22;
        h_2 = 3e-24;

    case 'OCXO'
        h_0 = 2e-25;
        h_1 = 7e-25;
        h_2 = 6e-25;

    case 'RUBIDIUM'
        h_0 = 2e-22;
        h_1 = 4.5e-26;
        h_2 = 1e-30;

    case 'CESIUM'
        h_0 = 2e-22;
        h_1 = 5e-27;
        h_2 = 1.5e-33;

    otherwise
        warning('Valid Receiver Clock Not Specified')
        h_0 = 0;
        h_1 = 0;
        h_2 = 0;
end

clkVar(1,1) = ((h_0/2)*timestep + 2*h_1*timestep^2 + (2/3)*pi^2*h_2*timestep^3);
clkVar(1,2) = (h_1*timestep + pi^2*h_2*timestep^2);
clkVar(2,1) = clkVar(1,2);
clkVar(2,2) = ((h_0/(2*timestep)) + 4*h_1 + (8/3)*pi^2 * h_2*timestep);

Qd(13:14,13:14) = clkVar;
end
