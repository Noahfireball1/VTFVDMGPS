function [rdot, vdot, omega_dot, euler_rates]  = CalcXdot(m,Ic_B,F_B,M_B,X)
% Constants
e = 0.0818191910428; % eccentricity
a = 6378137.0; % equatorial radius [meters]
omega_ie = 7.292115e-5; % Earth's rotation rate [rad/s]

u = X(1);
v = X(2);
w = X(3);
p = X(4);
q = X(5);
r = X(6);
lat = X(7);
long = X(8);
alt = X(9);
phi = X(10);
theta = X(11);
psi = X(12);

omega_skew = formskewsym([p;q;r]);
C_omega = [1 tan(theta)*sin(phi) tan(theta)*cos(phi);...
    0 cos(phi) -sin(phi);...
    0 sin(phi)/cos(theta) cos(phi)/cos(theta)];
C_n_b = [1 0 0;0 cos(phi) sin(phi); 0 -sin(phi) cos(phi)]*...
    [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)]*...
    [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1];
C_b_n = C_n_b';

R_N = (a*(1-e^2))/(1-e^2*sin(long)^2)^(3/2); % [checked]
R_E = (a)/(1 - e^2*sin(long)^2)^(1/2); % [checked]

omega_ie_n = [omega_ie*cos(long);0;-omega_ie*sin(long)]; % [checked]
omega_ie_n_skew = [0 -omega_ie_n(3) omega_ie_n(2); omega_ie_n(3) 0 -omega_ie_n(1); -omega_ie_n(2) omega_ie_n(1) 0];




rdot = [u/(R_N + alt);v/((R_E + alt)*cos(long));-w]; % [lat;long;alt] (radians,meters) Position derivative from earth to body in the nav frame

omega_en_n = [rdot(1)*cos(long);-rdot(2);rdot(1)*sin(long)]; % [checked]
omega_en_n_skew = [0 -omega_en_n(3) omega_en_n(2); omega_en_n(3) 0 -omega_en_n(1); -omega_en_n(2) omega_en_n(1) 0];

omega_nb_b = [p;q;r] - C_n_b*(omega_ie_n + omega_en_n);

vdot = C_b_n*(F_B/m) - (2*omega_ie_n_skew + omega_en_n_skew)*[u;v;w]; % [Nv;Ev;Dv] (meters) Velocity derivative from earth to body in the nav frame
euler_rates = C_omega*omega_nb_b ; % [phi_dot;theta_dot;psi_dot] (radians) euler rates from body to nav
omega_dot = Ic_B^-1*(M_B - omega_skew*(Ic_B*[p;q;r])); % [omega_x_dot, omega_y_dot, omega_z_dot] from inertial to body in the body frame

test = [vdot;omega_dot;rdot;euler_rates]*1/400 + X;

end