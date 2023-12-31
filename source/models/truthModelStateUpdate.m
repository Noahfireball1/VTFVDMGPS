function states = truthModelStateUpdate(forces,moments,oldStates,Time_Step,variance,clkVar)
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
R_N = (a*(1-e^2))/(1-e^2*sin(lat)^2)^(3/2);
R_E = (a)/(1 - e^2*sin(lat)^2)^(1/2);

omega_ie_n_skew = [0 sin(lat) 0;-sin(lat) 0 -cos(lat); 0 cos(lat) 0]*omega_ie;
omega_ie_n = [omega_ie*cos(lat);0;-omega_ie*sin(lat)];

% Linear Velocities
D = diag([1/(R_N + alt);1/((R_E + alt)*cos(lat));-1]);
rdot = D*([u;v;w] + sqrt(variance(7:9)).*randn(3,1)); % [lat;long;alt] (radians,meters) Position derivative from earth to body in the nav frame

omega_en_n = [v/(R_E + alt); -u/(R_N + alt); -v*tan(lat)/(R_E + alt)]; % [checked]
omega_en_n_skew = [0 -omega_en_n(3) omega_en_n(2); omega_en_n(3) 0 -omega_en_n(1); -omega_en_n(2) omega_en_n(1) 0];

omega_nb_b = [p;q;r] - C_n_b*(omega_ie_n + omega_en_n);
% Linear Accelerations
vdot = (C_b_n*(forces/m) - (2*omega_ie_n_skew + omega_en_n_skew)*[u;v;w]) + (sqrt(variance(1:3)).*randn(3,1)); % [Nv;Ev;Dv] (meters) Velocity derivative from earth to body in the nav frame

% Angular Accelerations
omega_dot = Ic_B^-1*(moments - omega_skew*(Ic_B*[p;q;r])) + sqrt(variance(4:6)).*randn(3,1); % [omega_x_dot, omega_y_dot, omega_z_dot] from inertial to body in the body frame

% Euler Rates
euler_rates = C_omega*omega_nb_b ; % [phi_dot;theta_dot;psi_dot] (radians) euler rates from body to nav



% Integrate
xDot = [vdot;omega_dot;rdot;euler_rates;0;0];
states = xDot*Time_Step + oldStates;

clkNoise = sqrt(clkVar)*randn(2,1);
states(13) = oldStates(13) + oldStates(14)*Time_Step + clkNoise(1);
states(14) = oldStates(14) + clkNoise(2);

end
