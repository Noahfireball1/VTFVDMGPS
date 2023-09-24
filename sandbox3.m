%% Formatting
clc
clear
close all
format short g

%% Sandbox

u = 75;
v = 0;
w = 0;
p = 0;
q = 0;
r = 0;
x = 0;
y = 0;
z = -1000;
phi = 0;
theta = 4*pi/180;
psi = 2*pi;
clkBias = 0;
clkDrift = 0;
timeStep = 1/400;
xk0 = [u;v;w;p;q;r;x;y;z;phi;theta;psi;clkBias;clkDrift];

m = 1202.7;
MOI = [4044.7 1994.9 5639.2];
cg = [-0.426;0;0.2];

A = zeros(14);
B = zeros(14,6);

f_b = [600 0 -1200];
m_b = f_b*10;

udot = f_b(1)/m + r*v - q*w;
vdot = f_b(2)/m + p*w - r*u;
wdot = f_b(3)/m + q*u - p*v;

pdot = (m_b(1) + MOI(2)*q*r - MOI(3)*q*r)/MOI(1);
qdot = (m_b(2) - MOI(1)*p*r + MOI(3)*p*r)/MOI(2);
rdot = (m_b(3) + MOI(1)*p*q - MOI(2)*p*q)/MOI(3);

xk0_dot = [udot;vdot;wdot;pdot;qdot;rdot];
A(1:3,1:3) = eye(3);
A(4:6,4:6) = eye(3);
A(1:3,7:9) = eye(3)*timeStep;
A(7:9,7:9) = eye(3);
A(4:6,10:12) = eye(3)*timeStep;
A(10:12,10:12) = eye(3);
A(13:14,13:14) = [1 timeStep; 0 1];

B(1:3,1:3) = timeStep;
B(4:6,4:6) = timeStep;
B(7:9,1:3) = timeStep^2/2;
B(10:12,4:6) = timeStep^2/2;

xk1 = A*xk0 + B*xk0_dot;

Q = diag(mean(B*xk0_dot*xk0_dot'*B'))