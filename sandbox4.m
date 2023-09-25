%% Formatting
clc
clear
close all
format short g

%% Sandbox

syms u v w p q r x y z phi theta ps Ixx Iyy Izz fb1 fb2 fb3 mb1 mb2 mb3 mass t

R_b_n = [1 0 0; 0 cos(phi) sin(phi);0 -sin(phi) cos(phi)]...
    *[cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)]...
    *[cos(ps) sin(ps) 0;-sin(ps) cos(ps) 0; 0 0 1];

R_omega = [1 tan(theta)*sin(phi) tan(theta)*cos(phi);0 cos(phi) -sin(phi);0 sin(phi)/cos(theta) cos(phi)/cos(theta)];

pdot = R_b_n*[u;v;w] + [x;y;z];
vdot = R_b_n*-[q*w - r*v;r*u - p*w;p*v - q*u];
omega = R_omega*[p;q;r];
omega_dot = [Iyy*q*r/Ixx - Izz*q*r/Ixx;...
Ixx*p*r/Iyy + Izz*p*r/Iyy;...
Ixx*p*q/Izz - Iyy*p*q/Izz];

f_n = R_b_n*[fb1;fb2;fb3];
m_n = R_b_n*[mb1;mb2;mb3];

w = [f_n./mass;m_n./[Ixx;Iyy;Izz]];

G = zeros(12,6);
G(1,1) = 1;
G(2,2) = 1;
G(3,3) = 1;
G(4,4) = 1;
G(5,5) = 1;
G(6,6) = 1;

Q = G*(w*w')*G'.*t;

% jacobian(F,[u v w p q r x y z phi theta ps])

%% du
dudot_du = -q*sin(theta) - r*cos(theta)*sin(ps);
dudot_dv = p*sin(theta) + r*cos(ps)*cos(theta);
dudot_dw = p*cos(theta)*sin(ps) - q*cos(ps)*cos(theta);
dudot_dp = v*sin(theta) + w*cos(theta)*sin(ps);
dudot_dq = -u*sin(theta) - w*cos(ps)*cos(theta);
dudot_dr = v*cos(ps)*cos(theta) - u*cos(theta)*sin(ps);
dudot_dx = 0;
dudot_dy = 0;
dudot_dz = 0;
dudot_dphi = 0;
dudot_dtheta = cos(theta)*(p*v - q*u) + cos(ps)*sin(theta)*(q*w - r*v) - sin(ps)*sin(theta)*(p*w - r*u);
dudot_dpsi = cos(theta)*sin(ps)*(q*w - r*v) + cos(ps)*cos(theta)*(p*w - r*u);
dudot_dclkBias = 0;
dudot_dclkDrift = 0;

%% dv
dvdot_du = q*cos(theta)*sin(phi) - r*(cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta));
dvdot_dv = - r*(cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta)) - p*cos(theta)*sin(phi);
dvdot_dw =  p*(cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta)) + q*(cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta));
dvdot_dp = w*(cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta)) - v*cos(theta)*sin(phi);
dvdot_dq = w*(cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta)) + u*cos(theta)*sin(phi);
dvdot_dr = - u*(cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta)) - v*(cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta));
dvdot_dx = 0;
dvdot_dy = 0;
dvdot_dz = 0;
dvdot_dphi = - (cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta))*(p*w - r*u) - (sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta))*(q*w - r*v) - cos(phi)*cos(theta)*(p*v - q*u);
dvdot_dtheta = sin(phi)*sin(theta)*(p*v - q*u) - cos(ps)*cos(theta)*sin(phi)*(q*w - r*v) + cos(theta)*sin(phi)*sin(ps)*(p*w - r*u);
dvdot_dpsi = (cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta))*(q*w - r*v) - (cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta))*(p*w - r*u);
dvdot_dclkBias = 0;
dvdot_dclkDrift = 0;

%% dw
dwdot_du = r*(cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta)) + q*cos(phi)*cos(theta);
dwdot_dv = r*(sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta)) - p*cos(phi)*cos(theta);
dwdot_dw = - p*(cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta)) - q*(sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta));
dwdot_dp = - w*(cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta)) - v*cos(phi)*cos(theta);
dwdot_dq = u*cos(phi)*cos(theta) - w*(sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta));
dwdot_dr = u*(cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta)) + v*(sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta));
dwdot_dx = 0;
dwdot_dy = 0;
dwdot_dz = 0;
dwdot_dphi = cos(theta)*sin(phi)*(p*v - q*u) - (cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta))*(q*w - r*v) - (cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta))*(p*w - r*u);
dwdot_dtheta = cos(phi)*sin(theta)*(p*v - q*u) - cos(phi)*cos(ps)*cos(theta)*(q*w - r*v) + cos(phi)*cos(theta)*sin(ps)*(p*w - r*u);
dwdot_dpsi = (sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta))*(p*w - r*u) - (cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta))*(q*w - r*v);
dwdot_dclkBias = 0;
dwdot_dclkDrift = 0;

%% dp
dpdot_du = 0;
dpdot_dv = 0;
dpdot_dw = 0;
dpdot_dp = 0;
dpdot_dq = (Iyy*r)/Ixx - (Izz*r)/Ixx;
dpdot_dr = (Iyy*q)/Ixx - (Izz*q)/Ixx;
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
dqdot_dp = (Ixx*r)/Iyy + (Izz*r)/Iyy;
dqdot_dq = 0;
dqdot_dr = (Ixx*p)/Iyy + (Izz*p)/Iyy;
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
drdot_dp = (Ixx*q)/Izz - (Iyy*q)/Izz;
drdot_dq = (Ixx*p)/Izz - (Iyy*p)/Izz;
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
dxdot_du = cos(ps)*cos(theta);
dxdot_dv = cos(theta)*sin(ps);
dxdot_dw = -sin(theta);
dxdot_dp = 0;
dxdot_dq = 0;
dxdot_dr = 0;
dxdot_dx = 0;
dxdot_dy = 0;
dxdot_dz = 0;
dxdot_dphi = 0;
dxdot_dtheta = - w*cos(theta) - u*cos(ps)*sin(theta) - v*sin(ps)*sin(theta);
dxdot_dpsi = v*cos(ps)*cos(theta) - u*cos(theta)*sin(ps);
dxdot_dclkBias = 0;
dxdot_dclkDrift = 0;

%% dy
dydot_du = cos(ps)*sin(phi)*sin(theta) - cos(phi)*sin(ps);
dydot_dv = cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta);
dydot_dw = cos(theta)*sin(phi);
dydot_dp = 0;
dydot_dq = 0;
dydot_dr = 0;
dydot_dx = 0;
dydot_dy = 0;
dydot_dz = 0;
dydot_dphi = u*(sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta)) - v*(cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta)) + w*cos(phi)*cos(theta);
dydot_dtheta = u*cos(ps)*cos(theta)*sin(phi) - w*sin(phi)*sin(theta) + v*cos(theta)*sin(phi)*sin(ps);
dydot_dpsi = - u*(cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta)) - v*(cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta));
dydot_dclkBias = 0;
dydot_dclkDrift = 0;

%% dz
dzdot_du = sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta);
dzdot_dv = cos(phi)*sin(ps)*sin(theta) - cos(ps)*sin(phi);
dzdot_dw = cos(phi)*cos(theta);
dzdot_dp = 0;
dzdot_dq = 0;
dzdot_dr = 0;
dzdot_dx = 0;
dzdot_dy = 0;
dzdot_dz = 0;
dzdot_dphi = u*(cos(phi)*sin(ps) - cos(ps)*sin(phi)*sin(theta)) - v*(cos(phi)*cos(ps) + sin(phi)*sin(ps)*sin(theta)) - w*cos(theta)*sin(phi);
dzdot_dtheta = u*cos(phi)*cos(ps)*cos(theta) - w*cos(phi)*sin(theta) + v*cos(phi)*cos(theta)*sin(ps);
dzdot_dpsi = u*(cos(ps)*sin(phi) - cos(phi)*sin(ps)*sin(theta)) + v*(sin(phi)*sin(ps) + cos(phi)*cos(ps)*sin(theta));
dzdot_dclkBias = 0;
dzdot_dclkDrift = 0;

%% dphi
dphidot_du = 0;
dphidot_dv = 0;
dphidot_dw = 0;
dphidot_dp = 1;
dphidot_dq = sin(phi)*tan(theta);
dphidot_dr = cos(phi)*tan(theta);
dphidot_dx = 0;
dphidot_dy = 0;
dphidot_dz = 0;
dphidot_dphi = q*cos(phi)*tan(theta) - r*sin(phi)*tan(theta);
dphidot_dtheta = r*cos(phi)*(tan(theta)^2 + 1) + q*sin(phi)*(tan(theta)^2 + 1);
dphidot_dpsi = 0;
dphidot_dclkBias = 0;
dphidot_dclkDrift = 0;

%% dtheta
dthetadot_du = 0;
dthetadot_dv = 0;
dthetadot_dw = 0;
dthetadot_dp = 0;
dthetadot_dq = cos(phi);
dthetadot_dr = -sin(phi);
dthetadot_dx = 0;
dthetadot_dy = 0;
dthetadot_dz = 0;
dthetadot_dphi = - r*cos(phi) - q*sin(phi);
dthetadot_dtheta = 0;
dthetadot_dpsi = 0;
dthetadot_dclkBias = 0;
dthetadot_dclkDrift = 0;

%% dpsi
dpsidot_du = 0;
dpsidot_dv = 0;
dpsidot_dw = 0;
dpsidot_dp = 0;
dpsidot_dq = sin(phi)/cos(theta);
dpsidot_dr = cos(phi)/cos(theta);
dpsidot_dx = 0;
dpsidot_dy = 0;
dpsidot_dz = 0;
dpsidot_dphi = (q*cos(phi))/cos(theta) - (r*sin(phi))/cos(theta);
dpsidot_dtheta = (r*cos(phi)*sin(theta))/cos(theta)^2 + (q*sin(phi)*sin(theta))/cos(theta)^2;
dpsidot_dpsi = 0;
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

