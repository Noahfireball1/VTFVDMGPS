%% Formatting
clc
clear
close all
format short g

%% Sandbox
e = 0.0818191910428; % eccentricity
a = 6378137.0; % equatorial radius [meters]
omega_ie = 7.292115e-5; % Earth's rotation rate [rad/s]

syms u v w p q r lat long alt phi theta ps Ixx Iyy Izz fb1 fb2 fb3 mb1 mb2 mb3 m omega_ie R_0 e clkBias clkDrift

% Constants


% u = X(1);
% v = X(2);
% w = X(3);
% p = X(4);
% q = X(5);
% r = X(6);
% lat = X(7);
% long = X(8);
% alt = X(9);
% phi = X(10);
% theta = X(11);
% ps = X(12);
forces = [fb1;fb2;fb3];

Ic_B = diag([Ixx Iyy Izz]);

moments = [mb1;mb2;mb3];
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
clkRates = [clkDrift;0];

F = [vdot;omega_dot;rdot;euler_rates];

J = simplify(jacobian(F,[u v w p q r lat long alt phi theta ps clkBias clkDrift]));
G = jacobian(F,[fb1 fb2 fb3 mb1 mb2 mb3]);


%% du
dudot_du = (v*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2));
dudot_dv = w/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2))) - 2*omega_ie*sin(long) + (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2));
dudot_dw = v/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)));
dudot_dp = 0;
dudot_dq = 0;
dudot_dr = 0;
dudot_dx = 0;
dudot_dy = v*((u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)) - 2*omega_ie*cos(long) + (3*a*e^2*u*cos(long)*sin(long)^2*(e^2 - 1))/((1 - e^2*sin(long)^2)^(5/2)*(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2)) + (v*w*sin(long))/(cos(long)^2*(alt + a/(1 - e^2*sin(long)^2)^(1/2))) - (a*e^2*v*w*sin(long))/((alt + a/(1 - e^2*sin(long)^2)^(1/2))^2*(1 - e^2*sin(long)^2)^(3/2));
dudot_dz = -(u*v*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2 - (v*w)/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2))^2);
dudot_dphi = 0;
dudot_dtheta = 0;
dudot_dpsi = 0;
dudot_dclkBias = 0;
dudot_dclkDrift = 0;

%% dv
dvdot_du = 2*omega_ie*sin(long) + (w*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)) - (2*u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2));
dvdot_dv = 0;
dvdot_dw =  2*omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2));
dvdot_dp = 0;
dvdot_dq = 0;
dvdot_dr = 0;
dvdot_dx = 0;
dvdot_dy = -u*((u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)) - 2*omega_ie*cos(long) + (3*a*e^2*u*cos(long)*sin(long)^2*(e^2 - 1))/((1 - e^2*sin(long)^2)^(5/2)*(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2)) - w*(2*omega_ie*sin(long) + (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)) - (3*a*e^2*u*cos(long)^2*sin(long)*(e^2 - 1))/((1 - e^2*sin(long)^2)^(5/2)*(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2));
dvdot_dz = (u^2*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2 - (u*w*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2;
dvdot_dphi = 0;
dvdot_dtheta = 0;
dvdot_dpsi = 0;
dvdot_dclkBias = 0;
dvdot_dclkDrift = 0;

%% dw
dwdot_du = -v/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2))) - (v*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2));
dwdot_dv = -2*omega_ie*cos(long) - u/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2))) - (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2));
dwdot_dw = 0;
dwdot_dp = 0;
dwdot_dq = 0;
dwdot_dr = 0;
dwdot_dx = 0;
dwdot_dy = v*(2*omega_ie*sin(long) + (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)) - (3*a*e^2*u*cos(long)^2*sin(long)*(e^2 - 1))/((1 - e^2*sin(long)^2)^(5/2)*(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2)) - (u*v*sin(long))/(cos(long)^2*(alt + a/(1 - e^2*sin(long)^2)^(1/2))) + (a*e^2*u*v*sin(long))/((alt + a/(1 - e^2*sin(long)^2)^(1/2))^2*(1 - e^2*sin(long)^2)^(3/2));
dwdot_dz = (u*v*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2 + (u*v)/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2))^2);
dwdot_dphi = 0;
dwdot_dtheta = 0;
dwdot_dpsi = 0;
dwdot_dclkBias = 0;
dwdot_dclkDrift = 0;

%% dp
dpdot_du = 0;
dpdot_dv = 0;
dpdot_dw = 0;
dpdot_dp = 0;
dpdot_dq = (Iyy*r - Izz*r)/Ixx;
dpdot_dr = (Iyy*q - Izz*q)/Ixx;
dpdot_dx = 0;
dpdot_dy = 0;
dpdot_dz = 0;
dpdot_dphi = 0;
dpdot_dtheta = 0;
dpdot_dpsi = 0;
dpdot_dclkBias = 0;
dpdot_dclkDrift = 0;

%% dq
dqdot_du = 0;
dqdot_dv = 0;
dqdot_dw = 0;
dqdot_dp = -(Ixx*r - Izz*r)/Iyy;
dqdot_dq = 0;
dqdot_dr = -(Ixx*p - Izz*p)/Iyy;
dqdot_dx = 0;
dqdot_dy = 0;
dqdot_dz = 0;
dqdot_dphi = 0;
dqdot_dtheta = 0;
dqdot_dpsi = 0;
dqdot_dclkBias = 0;
dqdot_dclkDrift = 0;

%% dr
drdot_du = 0;
drdot_dv = 0;
drdot_dw = 0;
drdot_dp = (Ixx*q - Iyy*q)/Izz;
drdot_dq = (Ixx*p - Iyy*p)/Izz;
drdot_dr = 0;
drdot_dx = 0;
drdot_dy = 0;
drdot_dz = 0;
drdot_dphi = 0;
drdot_dtheta = 0;
drdot_dpsi = 0;
drdot_dclkBias = 0;
drdot_dclkDrift = 0;

%% dx
dxdot_du = 1/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2));
dxdot_dv = 0;
dxdot_dw = 0;
dxdot_dp = 0;
dxdot_dq = 0;
dxdot_dr = 0;
dxdot_dx = 0;
dxdot_dy = (3*a*e^2*u*cos(long)*sin(long)*(e^2 - 1))/((1 - e^2*sin(long)^2)^(5/2)*(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2);
dxdot_dz = -u/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2;
dxdot_dphi = 0;
dxdot_dtheta = 0;
dxdot_dpsi = 0;
dxdot_dclkBias = 0;
dxdot_dclkDrift = 0;

%% dy
dydot_du = 0;
dydot_dv = 1/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)));
dydot_dw = 0;
dydot_dp = 0;
dydot_dq = 0;
dydot_dr = 0;
dydot_dx = 0;
dydot_dy = (v*sin(long))/(cos(long)^2*(alt + a/(1 - e^2*sin(long)^2)^(1/2))) - (a*e^2*v*sin(long))/((alt + a/(1 - e^2*sin(long)^2)^(1/2))^2*(1 - e^2*sin(long)^2)^(3/2));
dydot_dz = -v/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2))^2);
dydot_dphi = 0;
dydot_dtheta = 0;
dydot_dpsi = 0;
dydot_dclkBias = 0;
dydot_dclkDrift = 0;

%% dz
dzdot_du = 0;
dzdot_dv = 0;
dzdot_dw = -1;
dzdot_dp = 0;
dzdot_dq = 0;
dzdot_dr = 0;
dzdot_dx = 0;
dzdot_dy = 0;
dzdot_dz = 0;
dzdot_dphi = 0;
dzdot_dtheta = 0;
dzdot_dpsi = 0;
dzdot_dclkBias = 0;
dzdot_dclkDrift = 0;

%% dphi
dphidot_du = (sin(long)*sin(theta))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)) - cos(phi)*tan(theta)*((cos(long)*(sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta)))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)) + (cos(phi)*cos(theta)*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + sin(phi)*tan(theta)*((cos(long)*(cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta)))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)) - (cos(theta)*sin(long)*sin(phi))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) - (cos(long)*cos(ps)*cos(theta))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2));
dphidot_dv = (cos(theta)*sin(ps))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2))) - (cos(phi)*tan(theta)*(cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2))) + (sin(phi)*tan(theta)*(cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)));
dphidot_dw = 0;
dphidot_dp = 1;
dphidot_dq = sin(phi)*tan(theta);
dphidot_dr = cos(phi)*tan(theta);
dphidot_dx = 0;
dphidot_dy = sin(theta)*((u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)) - omega_ie*cos(long) + (3*a*e^2*u*cos(long)*sin(long)^2*(e^2 - 1))/((1 - e^2*sin(long)^2)^(5/2)*(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2)) + cos(ps)*cos(theta)*(omega_ie*sin(long) + (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)) - (3*a*e^2*u*cos(long)^2*sin(long)*(e^2 - 1))/((1 - e^2*sin(long)^2)^(5/2)*(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2)) + cos(phi)*tan(theta)*((sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta))*(omega_ie*sin(long) + (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)) - (3*a*e^2*u*cos(long)^2*sin(long)*(e^2 - 1))/((1 - e^2*sin(long)^2)^(5/2)*(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2)) - cos(phi)*cos(theta)*((u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)) - omega_ie*cos(long) + (3*a*e^2*u*cos(long)*sin(long)^2*(e^2 - 1))/((1 - e^2*sin(long)^2)^(5/2)*(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2)) - (v*sin(long)*(cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta)))/(cos(long)^2*(alt + a/(1 - e^2*sin(long)^2)^(1/2))) + (a*e^2*v*sin(long)*(cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta)))/((alt + a/(1 - e^2*sin(long)^2)^(1/2))^2*(1 - e^2*sin(long)^2)^(3/2))) - sin(phi)*tan(theta)*((cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta))*(omega_ie*sin(long) + (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)) - (3*a*e^2*u*cos(long)^2*sin(long)*(e^2 - 1))/((1 - e^2*sin(long)^2)^(5/2)*(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2)) + cos(theta)*sin(phi)*((u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)) - omega_ie*cos(long) + (3*a*e^2*u*cos(long)*sin(long)^2*(e^2 - 1))/((1 - e^2*sin(long)^2)^(5/2)*(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2)) - (v*sin(long)*(cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta)))/(cos(long)^2*(alt + a/(1 - e^2*sin(long)^2)^(1/2))) + (a*e^2*v*sin(long)*(cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta)))/((alt + a/(1 - e^2*sin(long)^2)^(1/2))^2*(1 - e^2*sin(long)^2)^(3/2))) + (v*cos(theta)*sin(long)*sin(ps))/(cos(long)^2*(alt + a/(1 - e^2*sin(long)^2)^(1/2))) - (a*e^2*v*cos(theta)*sin(long)*sin(ps))/((alt + a/(1 - e^2*sin(long)^2)^(1/2))^2*(1 - e^2*sin(long)^2)^(3/2));
dphidot_dz = cos(phi)*tan(theta)*((u*cos(long)*(sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta)))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2 + (v*(cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2))^2) + (u*cos(phi)*cos(theta)*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2) - sin(phi)*tan(theta)*((u*cos(long)*(cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta)))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2 + (v*(cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2))^2) - (u*cos(theta)*sin(long)*sin(phi))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2) - (u*sin(long)*sin(theta))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2 + (u*cos(long)*cos(ps)*cos(theta))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2 - (v*cos(theta)*sin(ps))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2))^2);
dphidot_dphi = cos(phi)*tan(theta)*(q + (cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta))*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + cos(theta)*sin(phi)*(omega_ie*sin(long) - (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + (v*(cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)))) - sin(phi)*tan(theta)*(r - (sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta))*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + cos(phi)*cos(theta)*(omega_ie*sin(long) - (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) - (v*(cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)))) - cos(phi)*tan(theta)*((cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta))*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + cos(theta)*sin(phi)*(omega_ie*sin(long) - (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + (v*(cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)))) - sin(phi)*tan(theta)*((sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta))*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) - cos(phi)*cos(theta)*(omega_ie*sin(long) - (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + (v*(cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2))));
dphidot_dtheta = cos(phi)*(tan(theta)^2 + 1)*(r - (sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta))*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + cos(phi)*cos(theta)*(omega_ie*sin(long) - (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) - (v*(cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)))) - cos(theta)*(omega_ie*sin(long) - (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) - cos(phi)*tan(theta)*(cos(phi)*sin(theta)*(omega_ie*sin(long) - (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + cos(phi)*cos(ps)*cos(theta)*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) - (v*cos(phi)*cos(theta)*sin(ps))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)))) + sin(phi)*(tan(theta)^2 + 1)*(q + (cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta))*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + cos(theta)*sin(phi)*(omega_ie*sin(long) - (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + (v*(cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)))) + cos(ps)*sin(theta)*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) - sin(phi)*tan(theta)*(sin(phi)*sin(theta)*(omega_ie*sin(long) - (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + cos(ps)*cos(theta)*sin(phi)*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) - (v*cos(theta)*sin(phi)*sin(ps))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)))) - (v*sin(ps)*sin(theta))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)));
dphidot_dpsi = sin(phi)*tan(theta)*((cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta))*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) - (v*(cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)))) - cos(phi)*tan(theta)*((cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta))*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) - (v*(sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)))) + cos(theta)*sin(ps)*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + (v*cos(ps)*cos(theta))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)));
dphidot_dclkBias = 0;
dphidot_dclkDrift = 0;

%% dtheta
dthetadot_du = sin(phi)*((cos(long)*(sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta)))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)) + (cos(phi)*cos(theta)*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + cos(phi)*((cos(long)*(cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta)))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)) - (cos(theta)*sin(long)*sin(phi))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)));
dthetadot_dv = (cos(phi)*(cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2))) + (sin(phi)*(cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)));
dthetadot_dw = 0;
dthetadot_dp = 0;
dthetadot_dq = cos(phi);
dthetadot_dr = -sin(phi);
dthetadot_dx = 0;
dthetadot_dy = - cos(phi)*((cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta))*(omega_ie*sin(long) + (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)) - (3*a*e^2*u*cos(long)^2*sin(long)*(e^2 - 1))/((1 - e^2*sin(long)^2)^(5/2)*(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2)) + cos(theta)*sin(phi)*((u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)) - omega_ie*cos(long) + (3*a*e^2*u*cos(long)*sin(long)^2*(e^2 - 1))/((1 - e^2*sin(long)^2)^(5/2)*(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2)) - (v*sin(long)*(cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta)))/(cos(long)^2*(alt + a/(1 - e^2*sin(long)^2)^(1/2))) + (a*e^2*v*sin(long)*(cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta)))/((alt + a/(1 - e^2*sin(long)^2)^(1/2))^2*(1 - e^2*sin(long)^2)^(3/2))) - sin(phi)*((sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta))*(omega_ie*sin(long) + (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)) - (3*a*e^2*u*cos(long)^2*sin(long)*(e^2 - 1))/((1 - e^2*sin(long)^2)^(5/2)*(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2)) - cos(phi)*cos(theta)*((u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)) - omega_ie*cos(long) + (3*a*e^2*u*cos(long)*sin(long)^2*(e^2 - 1))/((1 - e^2*sin(long)^2)^(5/2)*(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2)) - (v*sin(long)*(cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta)))/(cos(long)^2*(alt + a/(1 - e^2*sin(long)^2)^(1/2))) + (a*e^2*v*sin(long)*(cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta)))/((alt + a/(1 - e^2*sin(long)^2)^(1/2))^2*(1 - e^2*sin(long)^2)^(3/2)));
dthetadot_dz = - sin(phi)*((u*cos(long)*(sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta)))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2 + (v*(cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2))^2) + (u*cos(phi)*cos(theta)*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2) - cos(phi)*((u*cos(long)*(cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta)))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2 + (v*(cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2))^2) - (u*cos(theta)*sin(long)*sin(phi))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2);
dthetadot_dphi = sin(phi)*((cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta))*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + cos(theta)*sin(phi)*(omega_ie*sin(long) - (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + (v*(cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)))) - sin(phi)*(q + (cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta))*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + cos(theta)*sin(phi)*(omega_ie*sin(long) - (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + (v*(cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)))) - cos(phi)*((sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta))*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) - cos(phi)*cos(theta)*(omega_ie*sin(long) - (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + (v*(cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)))) - cos(phi)*(r - (sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta))*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + cos(phi)*cos(theta)*(omega_ie*sin(long) - (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) - (v*(cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2))));
dthetadot_dtheta = sin(phi)*(cos(phi)*sin(theta)*(omega_ie*sin(long) - (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + cos(phi)*cos(ps)*cos(theta)*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) - (v*cos(phi)*cos(theta)*sin(ps))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)))) - cos(phi)*(sin(phi)*sin(theta)*(omega_ie*sin(long) - (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + cos(ps)*cos(theta)*sin(phi)*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) - (v*cos(theta)*sin(phi)*sin(ps))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2))));
dthetadot_dpsi = cos(phi)*((cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta))*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) - (v*(cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)))) + sin(phi)*((cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta))*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) - (v*(sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2))));
dthetadot_dclkBias = 0;
dthetadot_dclkDrift = 0;

%% dpsi
dpsidot_du = (sin(phi)*((cos(long)*(cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta)))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)) - (cos(theta)*sin(long)*sin(phi))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))))/cos(theta) - (cos(phi)*((cos(long)*(sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta)))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)) + (cos(phi)*cos(theta)*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))))/cos(theta);
dpsidot_dv = (sin(phi)*(cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta)))/(cos(long)*cos(theta)*(alt + a/(1 - e^2*sin(long)^2)^(1/2))) - (cos(phi)*(cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta)))/(cos(long)*cos(theta)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)));
dpsidot_dw = 0;
dpsidot_dp = 0;
dpsidot_dq = sin(phi)/cos(theta);
dpsidot_dr = cos(phi)/cos(theta);
dpsidot_dx = 0;
dpsidot_dy = (cos(phi)*((sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta))*(omega_ie*sin(long) + (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)) - (3*a*e^2*u*cos(long)^2*sin(long)*(e^2 - 1))/((1 - e^2*sin(long)^2)^(5/2)*(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2)) - cos(phi)*cos(theta)*((u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)) - omega_ie*cos(long) + (3*a*e^2*u*cos(long)*sin(long)^2*(e^2 - 1))/((1 - e^2*sin(long)^2)^(5/2)*(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2)) - (v*sin(long)*(cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta)))/(cos(long)^2*(alt + a/(1 - e^2*sin(long)^2)^(1/2))) + (a*e^2*v*sin(long)*(cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta)))/((alt + a/(1 - e^2*sin(long)^2)^(1/2))^2*(1 - e^2*sin(long)^2)^(3/2))))/cos(theta) - (sin(phi)*((cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta))*(omega_ie*sin(long) + (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)) - (3*a*e^2*u*cos(long)^2*sin(long)*(e^2 - 1))/((1 - e^2*sin(long)^2)^(5/2)*(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2)) + cos(theta)*sin(phi)*((u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2)) - omega_ie*cos(long) + (3*a*e^2*u*cos(long)*sin(long)^2*(e^2 - 1))/((1 - e^2*sin(long)^2)^(5/2)*(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2)) - (v*sin(long)*(cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta)))/(cos(long)^2*(alt + a/(1 - e^2*sin(long)^2)^(1/2))) + (a*e^2*v*sin(long)*(cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta)))/((alt + a/(1 - e^2*sin(long)^2)^(1/2))^2*(1 - e^2*sin(long)^2)^(3/2))))/cos(theta);
dpsidot_dz = (cos(phi)*((u*cos(long)*(sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta)))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2 + (v*(cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2))^2) + (u*cos(phi)*cos(theta)*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2))/cos(theta) - (sin(phi)*((u*cos(long)*(cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta)))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2 + (v*(cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2))^2) - (u*cos(theta)*sin(long)*sin(phi))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))^2))/cos(theta);
dpsidot_dphi = (cos(phi)*(q + (cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta))*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + cos(theta)*sin(phi)*(omega_ie*sin(long) - (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + (v*(cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)))))/cos(theta) - (sin(phi)*(r - (sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta))*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + cos(phi)*cos(theta)*(omega_ie*sin(long) - (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) - (v*(cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)))))/cos(theta) - (cos(phi)*((cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta))*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + cos(theta)*sin(phi)*(omega_ie*sin(long) - (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + (v*(cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)))))/cos(theta) - (sin(phi)*((sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta))*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) - cos(phi)*cos(theta)*(omega_ie*sin(long) - (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + (v*(cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)))))/cos(theta);
dpsidot_dtheta = (cos(phi)*sin(theta)*(r - (sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta))*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + cos(phi)*cos(theta)*(omega_ie*sin(long) - (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) - (v*(cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)))))/cos(theta)^2 - (sin(phi)*(sin(phi)*sin(theta)*(omega_ie*sin(long) - (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + cos(ps)*cos(theta)*sin(phi)*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) - (v*cos(theta)*sin(phi)*sin(ps))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)))))/cos(theta) - (cos(phi)*(cos(phi)*sin(theta)*(omega_ie*sin(long) - (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + cos(phi)*cos(ps)*cos(theta)*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) - (v*cos(phi)*cos(theta)*sin(ps))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)))))/cos(theta) + (sin(phi)*sin(theta)*(q + (cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta))*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + cos(theta)*sin(phi)*(omega_ie*sin(long) - (u*sin(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) + (v*(cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)))))/cos(theta)^2;
dpsidot_dpsi = (sin(phi)*((cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta))*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) - (v*(cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)))))/cos(theta) - (cos(phi)*((cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta))*(omega_ie*cos(long) + (u*cos(long))/(alt - (a*(e^2 - 1))/(1 - e^2*sin(long)^2)^(3/2))) - (v*(sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta)))/(cos(long)*(alt + a/(1 - e^2*sin(long)^2)^(1/2)))))/cos(theta);
dpsidot_dclkBias = 0;
dpsidot_dclkDrift = 0;

%% dclkBias
dclkBiasdot_du = 0;
dclkBiasdot_dv = 0;
dclkBiasdot_dw = 0;
dclkBiasdot_dp = 0;
dclkBiasdot_dq = 0;
dclkBiasdot_dr = 0;
dclkBiasdot_dx = 0;
dclkBiasdot_dy = 0;
dclkBiasdot_dz = 0;
dclkBiasdot_dphi = 0;
dclkBiasdot_dtheta = 0;
dclkBiasdot_dpsi = 0;
dclkBiasdot_dclkBias = 0;
dclkBiasdot_dclkDrift = 1;

%% dclkDrift
dclkDriftdot_du = 0;
dclkDriftdot_dv = 0;
dclkDriftdot_dw = 0;
dclkDriftdot_dp = 0;
dclkDriftdot_dq = 0;
dclkDriftdot_dr = 0;
dclkDriftdot_dx = 0;
dclkDriftdot_dy = 0;
dclkDriftdot_dz = 0;
dclkDriftdot_dphi = 0;
dclkDriftdot_dtheta = 0;
dclkDriftdot_dpsi = 0;
dclkDriftdot_dclkBias = 0;
dclkDriftdot_dclkDrift = 0;

F = [dudot_du dudot_dv dudot_dw dudot_dp dudot_dq dudot_dr dudot_dx dudot_dy dudot_dz dudot_dphi dudot_dtheta dudot_dpsi dudot_dclkBias dudot_dclkDrift;...
     dvdot_du dvdot_dv dvdot_dw dvdot_dp dvdot_dq dvdot_dr dvdot_dx dvdot_dy dvdot_dz dvdot_dphi dvdot_dtheta dvdot_dpsi dvdot_dclkBias dvdot_dclkDrift;...
     dwdot_du dwdot_dv dwdot_dw dwdot_dp dwdot_dq dwdot_dr dwdot_dx dwdot_dy dwdot_dz dwdot_dphi dwdot_dtheta dwdot_dpsi dwdot_dclkBias dwdot_dclkDrift;...
     dpdot_du dpdot_dv dpdot_dw dpdot_dp dpdot_dq dpdot_dr dpdot_dx dpdot_dy dpdot_dz dpdot_dphi dpdot_dtheta dpdot_dpsi dpdot_dclkBias dpdot_dclkDrift;...
     dqdot_du dqdot_dv dqdot_dw dqdot_dp dqdot_dq dqdot_dr dqdot_dx dqdot_dy dqdot_dz dqdot_dphi dqdot_dtheta dqdot_dpsi dqdot_dclkBias dqdot_dclkDrift;...
     drdot_du drdot_dv drdot_dw drdot_dp drdot_dq drdot_dr drdot_dx drdot_dy drdot_dz drdot_dphi drdot_dtheta drdot_dpsi drdot_dclkBias drdot_dclkDrift;...
     dxdot_du dxdot_dv dxdot_dw dxdot_dp dxdot_dq dxdot_dr dxdot_dx dxdot_dy dxdot_dz dxdot_dphi dxdot_dtheta dxdot_dpsi dxdot_dclkBias dxdot_dclkDrift;...
     dydot_du dydot_dv dydot_dw dydot_dp dydot_dq dydot_dr dydot_dx dydot_dy dydot_dz dydot_dphi dydot_dtheta dydot_dpsi dydot_dclkBias dydot_dclkDrift;...
     dzdot_du dzdot_dv dzdot_dw dzdot_dp dzdot_dq dzdot_dr dzdot_dx dzdot_dy dzdot_dz dzdot_dphi dzdot_dtheta dzdot_dpsi dzdot_dclkBias dzdot_dclkDrift;...
     dphidot_du dphidot_dv dphidot_dw dphidot_dp dphidot_dq dphidot_dr dphidot_dx dphidot_dy dphidot_dz dphidot_dphi dphidot_dtheta dphidot_dpsi dphidot_dclkBias dphidot_dclkDrift;...
     dthetadot_du dthetadot_dv dthetadot_dw dthetadot_dp dthetadot_dq dthetadot_dr dthetadot_dx dthetadot_dy dthetadot_dz dthetadot_dphi dthetadot_dtheta dthetadot_dpsi dthetadot_dclkBias dthetadot_dclkDrift;...
     dpsidot_du dpsidot_dv dpsidot_dw dpsidot_dp dpsidot_dq dpsidot_dr dpsidot_dx dpsidot_dy dpsidot_dz dpsidot_dphi dpsidot_dtheta dpsidot_dpsi dpsidot_dclkBias dpsidot_dclkDrift;...
     dclkBiasdot_du dclkBiasdot_dv dclkBiasdot_dw dclkBiasdot_dp dclkBiasdot_dq dclkBiasdot_dr dclkBiasdot_dx dclkBiasdot_dy dclkBiasdot_dz dclkBiasdot_dphi dclkBiasdot_dtheta dclkBiasdot_dpsi dclkBiasdot_dclkBias dclkBiasdot_dclkDrift;...
     dclkDriftdot_du dclkDriftdot_dv dclkDriftdot_dw dclkDriftdot_dp dclkDriftdot_dq dclkDriftdot_dr dclkDriftdot_dx dclkDriftdot_dy dclkDriftdot_dz dclkDriftdot_dphi dclkDriftdot_dtheta dclkDriftdot_dpsi dclkDriftdot_dclkBias dclkDriftdot_dclkDrift];
J-F
