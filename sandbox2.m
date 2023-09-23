%% Formatting
clc
clear
close all
format short g

%% Sandbox
orient = [30 20 135]*pi/180;
refLLA = [-34.9 138.5 10000];
lat = refLLA(1);
long = refLLA(2);
alt = refLLA(3);

% 1. Convert LLA Position to ECEF
ECEF = lla2ecef(refLLA);

% 2. Find Local NED Axes
N0 = [0;0;1];
E0 = [0;1;0];
D0 = [-1;0;0];

E = whatisR(refLLA(2)*pi/180,N0)*E0;
N = whatisR(refLLA(1)*pi/180,-E)*N0;
D = cross(N,E);

X0 = N;
Y0 = E;
Z0 = D;

X1 = whatisR(orient(3),Z0)*X0;
Y1 = whatisR(orient(3),Z0)*Y0;
Z1 = whatisR(orient(3),Z0)*Z0;

X2 = whatisR(orient(2),Y1)*X1;
Y2 = whatisR(orient(2),Y1)*Y1;
Z2 = whatisR(orient(2),Y1)*Z1;

X3 = whatisR(orient(1),X2)*X2
Y3 = whatisR(orient(1),X2)*Y2
Z3 = whatisR(orient(1),X2)*Z2


% 4. Calculate Euler Angles in New Frame
x0 = [1;0;0];
y0 = [0;1;0];
z0 = [0;0;1];

psi = (atan2( dot(X3,y0),dot(X3,x0) ));
theta = (atan2(-dot(X3,z0), sqrt((dot(X3,x0))^2+(dot(X3,y0))^2) ));


y2 = whatisR(psi,z0)*y0;
z2 = whatisR(theta,y2)*z0;

psi*180/pi
theta*180/pi
phi = (atan2( dot(Y3,z2), dot(Y3,y2) ))*(180/pi)

function Rn0 = whatisR(theta,N)
n1 = N(1);
n2 = N(2);
n3 = N(3);

R = N*N';
R1 = [0 -n3 n2;n3 0 -n1;-n2 n1 0];
Rn0 = (1 - cos(theta))*R + cos(theta)*eye(3) + sin(theta)*R1;

end

