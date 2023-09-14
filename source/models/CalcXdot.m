function [X_Dot,predictedCovariance]  = CalcXdot(m,Ic_B,F_B,M_B,X,CG,InvMassMatrix,timestep,oldCovariance,Q)

% muVector = zeros(3,1);

Xcg = CG(1);
Ycg = CG(2);
Zcg = CG(3);

CGssm = [0 -Zcg Ycg;Zcg 0 -Xcg;-Ycg Xcg 0];

syms u v w p q r phi theta Xpsi

lambdaX = F_B(1)/m + r*v - q*w - 0 - 2*q*0 + 2*r*0 + Xcg*(q^2 + r^2) - Ycg*p*q - Zcg*p*r;
lambdaY = F_B(2)/m + p*w - r*u - 0 - 2*r*0 + 2*p*0 - Xcg*p*q + Ycg*(p^2 +r^2) - Zcg*q*r;
lambdaZ = F_B(3)/m + q*u - p*v - 0 - 2*p*0 + 2*q*0 - Xcg*p*r - Ycg*q*r + Zcg*(p^2 + q^2);

MomentVector = [p;q;r];

MomentVectorssm = [0 -r q;r 0 -p;-q p 0];

muVector = [M_B(1);M_B(2);M_B(3)] - (MomentVectorssm*Ic_B*MomentVector) - m*CGssm*MomentVectorssm*[u;v;w];
LAMEparams = [lambdaX;lambdaY;lambdaZ;muVector];
X_Dot(1,1:6) = InvMassMatrix*LAMEparams;

xdot = cos(theta)*cos(Xpsi)*u + sin(phi)*sin(theta)*cos(Xpsi) - cos(phi)*sin(Xpsi)*v + cos(phi)*sin(theta)*cos(Xpsi) + sin(phi)*sin(Xpsi)*w;
ydot = cos(theta)*sin(Xpsi)*u + sin(phi)*sin(theta)*sin(Xpsi) + cos(phi)*cos(Xpsi)*v + cos(phi)*sin(theta)*sin(Xpsi) - sin(phi)*cos(Xpsi)*w;
zdot = -sin(theta)*u + sin(phi)*cos(theta)*v + cos(phi)*cos(theta)*w;
phidot = p + q*sin(phi)*tan(theta) + r*cos(phi)*tan(theta);
thdot = q*cos(phi) - r*sin(phi);
psidot = q*sin(phi)/cos(theta) + r*cos(phi)/cos(theta);
X_Dot(1,7:12) = [xdot;ydot;zdot;phidot;thdot;psidot];

bigF = jacobian(X_Dot);

bigF = subs(bigF,{u, v, w, p, q, r, phi, theta, Xpsi},{X(1), X(2), X(3), X(4), X(5), X(6), X(10), X(11), X(12)});

finalF = double([bigF(1,1) bigF(1,2) bigF(1,3) bigF(1,4) bigF(1,5) bigF(1,6) 0 0 0 bigF(1,7) bigF(1,8) bigF(1,9);...
          bigF(2,1) bigF(2,2) bigF(2,3) bigF(2,4) bigF(2,5) bigF(2,6) 0 0 0 bigF(2,7) bigF(2,8) bigF(2,9);...
          bigF(3,1) bigF(3,2) bigF(3,3) bigF(3,4) bigF(3,5) bigF(3,6) 0 0 0 bigF(3,7) bigF(3,8) bigF(3,9);...
          bigF(4,1) bigF(4,2) bigF(4,3) bigF(4,4) bigF(4,5) bigF(4,6) 0 0 0 bigF(4,7) bigF(4,8) bigF(4,9);...
          bigF(5,1) bigF(5,2) bigF(5,3) bigF(5,4) bigF(5,5) bigF(5,6) 0 0 0 bigF(5,7) bigF(5,8) bigF(5,9);...
          bigF(6,1) bigF(6,2) bigF(6,3) bigF(6,4) bigF(6,5) bigF(6,6) 0 0 0 bigF(6,7) bigF(6,8) bigF(6,9);...
          0         0         0         0         0         0         1 0 0 0          0          0         ;...
          0         0         0         0         0         0         0 1 0 0          0          0         ;...
          0         0         0         0         0         0         0 0 1 0          0          0         ;...
          bigF(10,1) bigF(10,2) bigF(10,3) bigF(10,4) bigF(10,5) bigF(10,6) 0 0 0 bigF(10,7) bigF(10,8) bigF(10,9);...
          bigF(11,1) bigF(11,2) bigF(11,3) bigF(11,4) bigF(11,5) bigF(11,6) 0 0 0 bigF(11,7) bigF(11,8) bigF(11,9);...
          bigF(12,1) bigF(12,2) bigF(12,3) bigF(12,4) bigF(12,5) bigF(12,6) 0 0 0 bigF(12,7) bigF(12,8) bigF(12,9)]);...


bigPHI = expm(finalF.*timestep);

predictedCovariance = bigPHI*oldCovariance*bigPHI' + Q;
X_Dot = double(subs(X_Dot,{u, v, w, p, q, r, phi, theta, Xpsi},{X(1), X(2), X(3), X(4), X(5), X(6), X(10), X(11), X(12)}));

stop = 1;
end

