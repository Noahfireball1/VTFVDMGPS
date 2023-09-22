%% Formatting
clc
clear
close all
format short g

%% Sandbox

syms m ...
    Ic_B ...
    fb1 fb2 fb3 ...
    mb1 mb2 mb3 ...
    Xcg Ycg Zcg ...
    imm11 imm12 imm13 imm14 imm15 imm16...
    imm21 imm22 imm23 imm24 imm25 imm26...
    imm31 imm32 imm33 imm34 imm35 imm36...
    imm41 imm42 imm43 imm44 imm45 imm46...
    imm51 imm52 imm53 imm54 imm55 imm56...
    imm61 imm62 imm63 imm64 imm65 imm66 ...
    u v w p q r x y z phi theta ps clkBias clkDrift

InvMassMatrix = [ imm11 imm12 imm13 imm14 imm15 imm16;...
    imm21 imm22 imm23 imm24 imm25 imm26;...
    imm31 imm32 imm33 imm34 imm35 imm36;...
    imm41 imm42 imm43 imm44 imm45 imm46;...
    imm51 imm52 imm53 imm54 imm55 imm56;...
    imm61 imm62 imm63 imm64 imm65 imm66];


CGssm = [0 -Zcg Ycg;Zcg 0 -Xcg;-Ycg Xcg 0];


lambdaX = fb1/m + r*v - q*w - 0 - 2*q*0 + 2*r*0 + Xcg*(q^2 + r^2) - Ycg*p*q - Zcg*p*r;
lambdaY = fb2/m + p*w - r*u - 0 - 2*r*0 + 2*p*0 - Xcg*p*q + Ycg*(p^2 +r^2) - Zcg*q*r;
lambdaZ = fb3/m + q*u - p*v - 0 - 2*p*0 + 2*q*0 - Xcg*p*r - Ycg*q*r + Zcg*(p^2 + q^2);

MomentVector = [p;q;r];

MomentVectorssm = [0 -r q;r 0 -p;-q p 0];

muVector = [mb1;mb2;mb3] - (MomentVectorssm*Ic_B*MomentVector) - m*CGssm*MomentVectorssm*[u;v;w];
LAMEparams = [lambdaX;lambdaY;lambdaZ;muVector];
X_Dot(1:6,:) = InvMassMatrix*LAMEparams;

cphi = cos(phi);
sphi = sin(phi);

ctheta = cos(theta);
stheta = sin(theta);
ttheta = tan(theta);

cpsi = cos(ps);
spsi = sin(ps);

xdot = ctheta*cpsi*u + sphi*stheta*cpsi - cphi*spsi*v + cphi*stheta*cpsi + sphi*spsi*w + x;
ydot = ctheta*spsi*u + sphi*stheta*spsi + cphi*cpsi*v + cphi*stheta*spsi - sphi*cpsi*w  + y;
zdot = -stheta*u + sphi*ctheta*v + cphi*ctheta*w + z;
phidot = p + q*sphi*ttheta + r*cphi*ttheta;
thdot = q*cphi - r*sphi;
psidot = q*sphi/ctheta + r*cphi/ctheta;
X_Dot(7:12,:) = [xdot;ydot;zdot;phidot;thdot;psidot];

X_Dot(13:14,:) = [clkDrift;0];

jacobian(X_Dot,[u,v,w,p,q,r,x,y,z,phi,theta,ps,clkBias,clkDrift])