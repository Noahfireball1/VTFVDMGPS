function X_Dot = CalcXdot(m,Ic_B,F_B,M_B,X,CG,InvMassMatrix,timeStep)

u = X(1);
v = X(2);
w = X(3);
p = X(4);
q = X(5);
r = X(6);
x = X(7);
y = X(8);
z = X(9);
phi = X(10);
theta = X(11);
psi = X(12);


Xcg = CG(1);
Ycg = CG(2);
Zcg = CG(3);

CGssm = [0 -Zcg Ycg;Zcg 0 -Xcg;-Ycg Xcg 0];


lambdaX = F_B(1)/m + r*v - q*w - 0 - 2*q*0 + 2*r*0 + Xcg*(q^2 + r^2) - Ycg*p*q - Zcg*p*r;
lambdaY = F_B(2)/m + p*w - r*u - 0 - 2*r*0 + 2*p*0 - Xcg*p*q + Ycg*(p^2 +r^2) - Zcg*q*r;
lambdaZ = F_B(3)/m + q*u - p*v - 0 - 2*p*0 + 2*q*0 - Xcg*p*r - Ycg*q*r + Zcg*(p^2 + q^2);

MomentVector = [p;q;r];

MomentVectorssm = [0 -r q;r 0 -p;-q p 0];

muVector = [M_B(1);M_B(2);M_B(3)] - (MomentVectorssm*Ic_B*MomentVector) - m*CGssm*MomentVectorssm*[u;v;w];
LAMEparams = [lambdaX;lambdaY;lambdaZ;muVector];
X_Dot(1:6,:) = InvMassMatrix*LAMEparams;

cphi = cos(phi);
sphi = sin(phi);

ctheta = cos(theta);
stheta = sin(theta);
ttheta = tan(theta);

cpsi = cos(psi);
spsi = sin(psi);

xdot = ctheta*cpsi*u + sphi*stheta*cpsi - cphi*spsi*v + cphi*stheta*cpsi + sphi*spsi*w;
ydot = ctheta*spsi*u + sphi*stheta*spsi + cphi*cpsi*v + cphi*stheta*spsi - sphi*cpsi*w;
zdot = -stheta*u + sphi*ctheta*v + cphi*ctheta*w;
phidot = p + q*sphi*ttheta + r*cphi*ttheta;
thdot = q*cphi - r*sphi;
psidot = q*sphi/ctheta + r*cphi/ctheta;
X_Dot(7:12,:) = [xdot;ydot;zdot;phidot;thdot;psidot];

X_Dot(13:14,:) = [X(14);0];
X_Dot = X_Dot + randn(14,1).*[0.15 0.15 0.15 0 0 0 1.5 1.5 1.5 0 0 0 1e-10 0]';

end