function states = equationsOfMotionWithNoise(oldStates,forces,moments,Time_Step,variance,clkVar)
% Constants
e = 0.0818191910428; % eccentricity
R_0 = 6378137.0; % equatorial radius [meters]
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

omega_skew = [0 -r q; r 0 -p; -q p 0];

% Rotation Matrix for Rotating Angular Rates into Euler Rates (Groves Eq. 2.58)
C_omega = [1 tan(theta)*sin(phi) tan(theta)*cos(phi);0 cos(phi) -sin(phi);0 sin(phi)/cos(theta) cos(phi)/cos(theta)];

% Rotation Matrix for Rotating Vectors from Navigation Frame to Body Frame (Groves Eq. 2.22)
C_n_b = [1 0 0;0 cos(phi) sin(phi); 0 -sin(phi) cos(phi)]*...
    [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)]*...
    [cos(ps) sin(ps) 0; -sin(ps) cos(ps) 0; 0 0 1];
C_b_n = C_n_b';


% Meridian Radius of Curvature (Groves Eq. 2.105)
R_N = (R_0*(1-e^2))/(1 - (e^2)*sin(lat)^2)^(3/2);

% Tranverse Radius of Curvature (Groves Eq. 2.106)
R_E = (R_0)/(sqrt(1 - (e^2)*sin(lat)^2));

% Earth's Rotation with Respect to the Navigation Frame in Skew Form (Groves Eq. 5.41)
omega_ie_n_skew = [0 sin(lat) 0;...
    -sin(lat) 0 -cos(lat);...
    0 cos(lat) 0]*omega_ie;

% Earth's Rotaion with Respect to the Navigation Frame in Vector Form (Groves 2.123)
omega_ie_n = [omega_ie*cos(lat);...
    0;...
    -omega_ie*sin(lat)];

% Transport Rate with Respect to the Navigation Frame in Vector Form (Groves Eq. 5.44)
omega_en_n = [v/(R_E + alt);...
    -u/(R_N + alt);...
    (-v*tan(lat))/(R_E + alt)];

% Transport Rate with Respect to the Navigation Frame in Skew Form (Groves Eq. 5.44)
omega_en_n_skew = [0 -omega_en_n(3) omega_en_n(2);...
    omega_en_n(3) 0 -omega_en_n(1);...
    -omega_en_n(2) omega_en_n(1) 0];

% Curvilinear Position Derivatives (Groves Eq. 2.111) [radians]
rdot = [u/(R_N + alt);...
    v/((R_E + alt)*cos(lat));...
    -w];

% Euler Angle Derivatives (Groves Eq. 5.45) [radians]
omega_nb_b = [p;q;r] - C_n_b*(omega_ie_n + omega_en_n);

% Velocity Derivates (Groves Eq. 5.53)
vdot = -(omega_en_n_skew + 2*omega_ie_n_skew)*[u;v;w] + C_b_n*(forces/m);

% Angular Rate Derivatives with Respect to the Body Frame (Khaghani 2016, Eq. 29)
omega_dot = Ic_B^-1*(moments - omega_skew*(Ic_B*[p;q;r]));

% Euler Angle Derivatives with Respect to the Navigation Frame (Groves Eq. 2.58)
euler_rates = C_omega*omega_nb_b ;

% Clock Derivatives (Simple Math)
clkRates = [oldStates(14);0];

% Integrate and Add Noise
noiseTerms = blkdiag(diag(variance),clkVar);
xDot = [vdot;omega_dot;rdot;euler_rates;clkRates] + sqrt(noiseTerms)*randn(14,1);
states = xDot*Time_Step + oldStates;
end

