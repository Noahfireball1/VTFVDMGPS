function [T,Q,eta_p,dOmega_dt] = CTCQPropeff(Prop_Pitch,power,Q_drive,rho,V,n,PROP)

%% Propellor Constants
d = 1.9; %m
I_prop = 5; %kg*m^2

%% Governor State Vector Variable Extraction
beta = Prop_Pitch;
% beta_dot = Gov_States(2);
% P_oil = Gov_States(3);
% x = Gov_States(4);

%% Gridded Interpolant Setup
if isempty(PROP)
    global Prop;
else
    Prop = PROP;
end

%% Propellor Dynamics (Thrust and Torque)
J = (V*cosd(beta))/(n*d);

C_T = 0.0041*J^3 - 0.0353*J^2 - 0.0177*J + 0.0947; %thrust coefficient
C_Q = 0.001*J^4 - 0.008*J^3 + 0.0078*J^2 - 0.0035*J + 0.0114; %torque coefficient

T = C_T*rho*n^2*d^4;
Q = C_Q*rho*n^2*d^5;
denom = (T/(pi*d/4*V^2*rho/2)) + 1;
eta_p = (2)/(1 + sqrt(denom));

%----------RPS-ODE----------
dOmega_dt = (30/pi)*((Q_drive - Q)/I_prop);

end