function [states,Qd] = predictStates(oldStates,forces,moments,Time_Step,variance,Phi,clkVar)
% Constants
e = 0.0818191910428; % eccentricity
a = 6378137.0; % equatorial radius [meters]
omega_ie = 7.292115e-5; % Earth's rotation rate [rad/s]
m = 1202.7; % Mass of Vehicle [kg]
Ixx = 4044.7; % Mass Moments of Inertia
Iyy = 1994.9;
Izz = 5639.2;
Ic_B = [Ixx 0 0; 0 Iyy 0; 0 0 Izz];

% Defining Each State Index
u = oldStates(1);
v = oldStates(2);
w = oldStates(3);
p = oldStates(4);
q = oldStates(5);
r = oldStates(6);
lat = oldStates(7);
long = oldStates(8);
alt = oldStates(9);
phi = oldStates(10);
theta = oldStates(11);
ps = oldStates(12);

% Rotation Matrices
C_omega = [1 tan(theta)*sin(phi) tan(theta)*cos(phi);...
    0 cos(phi) -sin(phi);...
    0 sin(phi)/cos(theta) cos(phi)/cos(theta)];
C_n_b = [1 0 0;0 cos(phi) sin(phi); 0 -sin(phi) cos(phi)]*...
    [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)]*...
    [cos(ps) sin(ps) 0; -sin(ps) cos(ps) 0; 0 0 1];
C_b_n = C_n_b';

% Pre-processing Calculations
omega_skew = [0 -r q; r 0 -p; -q p 0];
R_N = (a*(1-e^2))/(1-e^2*sin(long)^2)^(3/2);
R_E = (a)/(1 - e^2*sin(long)^2)^(1/2);

omega_ie_n = [omega_ie*cos(long);0;-omega_ie*sin(long)]; % [checked]
omega_ie_n_skew = [0 -omega_ie_n(3) omega_ie_n(2); omega_ie_n(3) 0 -omega_ie_n(1); -omega_ie_n(2) omega_ie_n(1) 0];

% Linear Velocities
rDot = [u/(R_N + alt);v/((R_E + alt)*cos(long));-w]; % [lat;long;alt] (radians,meters) Position derivative from earth to body in the nav frame

omega_en_n = [rDot(1)*cos(long);-rDot(2);rDot(1)*sin(long)]; % [checked]
omega_en_n_skew = [0 -omega_en_n(3) omega_en_n(2); omega_en_n(3) 0 -omega_en_n(1); -omega_en_n(2) omega_en_n(1) 0];

omega_nb_b = [p;q;r] - C_n_b*(omega_ie_n + omega_en_n);
% Linear Accelerations
vDot = C_b_n*(forces/m) - (2*omega_ie_n_skew + omega_en_n_skew)*[u;v;w]; % [Nv;Ev;Dv] (meters) Velocity derivative from earth to body in the nav frame

% Angular Accelerations
omegaDot = Ic_B^-1*(moments - omega_skew*(Ic_B*[p;q;r])); % [omega_x_dot, omega_y_dot, omega_z_dot] from inertial to body in the body frame

% Euler Rates
eulerRates = C_omega*omega_nb_b ; % [phi_dot;theta_dot;psi_dot] (radians) euler rates from body to nav

% Clock term rates
clkNoise = sqrt(clkVar)*randn(2,1);
clkRates = [oldStates(14) + clkNoise(1);clkNoise(2)];

% Adding Noise (G*w)
G = zeros(12,6);
G(1:6,1:6) = eye(6);
w = sqrt(variance);

noiseXdot = [vDot;omegaDot;rDot;eulerRates]  + G*w.*randn(12,1);

Qd = Phi(1:12,1:12)*G*diag(variance)*G';
Qd = blkdiag(Qd,clkVar);
% Integrate
xDot = [noiseXdot;clkRates];
states = xDot*Time_Step + oldStates;
end

