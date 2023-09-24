function newCov = calculateCovariance(xk0,timeStep,oldCov,refLL)
newu = xk0(1);
newv = xk0(2);
neww = xk0(3);
newp = xk0(4);
newq = xk0(5);
newr = xk0(6);
newx = xk0(7);
newy = xk0(8);
newz = xk0(9);
newphi = xk0(10);
newtheta = xk0(11);
newps = xk0(12);
clkBias = xk0(13);
clkDrift = xk0(14);

cphi = cos(newphi);
sphi = sin(newphi);

stheta = sin(newtheta);
ctheta = cos(newtheta);

cpsi = cos(newps);
spsi = sin(newps);

m = 1202.7;
MOI = [4044.7 1994.9 5639.2];

LLA = flat2lla([newx newy newz],refLL,0,0,'WGS84');

C_l2e = [-cosd(LLA(2))*sind(LLA(1)) -sind(LLA(2))*sind(LLA(1)) cosd(LLA(1));
        sind(LLA(2)) cosd(LLA(2)) 0;
        cosd(LLA(2))*cosd(LLA(1))  cosd(LLA(1))*sind(LLA(2)) sind(LLA(1));
        ]';

R_1 = [1 0 0; 0 cphi sphi; 0 -sphi cphi];
R_2 = [ctheta 0 -sphi; 0 1 0; stheta 0 ctheta];
R_3 = [cpsi spsi 0; -spsi cpsi 0; 0 0 1];

C_b2l = R_1*R_2*R_3;

C_b2e = C_l2e*C_b2l;



syms u v w p q r x y z phi theta ps

acc_b = [r*v - q*w;...
    p*w - r*u;...
    q*u - p*v];

acc_e = C_b2e*acc_b;

omega_b = [(MOI(2)*q*r - MOI(3)*q*r)/MOI(1);...
(MOI(1)*p*r + MOI(3)*p*r)/MOI(2);...
(MOI(1)*p*q - MOI(2)*p*q)/MOI(3)];

omega_e = C_b2e*omega_b; 

DCM = [cos(theta)*cos(ps), sin(phi)*sin(theta)*cos(ps) - cos(phi)*sin(ps), cos(phi)+sin(theta)*cos(ps) + sin(phi)*sin(ps);...
    cos(theta)*sin(ps), sin(phi)*sin(theta)*sin(ps) + cos(phi)*cos(ps), cos(phi)*sin(theta)*sin(ps) - sin(phi)*cos(ps);...
    -sin(theta)    , sin(phi)*cos(theta)                 , cos(phi)*cos(theta)];

DCM1 = [1 sin(phi)*sin(theta)/cos(theta) cos(phi)*sin(theta)/cos(theta);...
    0 cos(phi) -sin(phi);...
    0 sin(phi)/cos(theta) cos(phi)/cos(theta)];

xyzdot = DCM*[u;v;w] + [x;y;z];

xyzdot_e = C_l2e*xyzdot;

ptpdot = DCM1*[p;q;r];

ptpdot_e = C_l2e*ptpdot;


xk0_dot = [acc_e;omega_e;xyzdot_e;ptpdot_e];

F = jacobian(xk0_dot,[u v w p q r x y z phi theta ps]);

F_actual = subs(F,[u v w p q r x y z phi theta ps]...
    ,[newu newv neww newp newq newr newx newy newz newphi newtheta newps]);

F_actual = double(vpa(F_actual)) ;

PHI = eye(12) + F_actual*timeStep;
PHI(13:14,13:14) = [1 timeStep; 0 1];
PHI(7:9,7:9) = PHI(7:9,7:9) - eye(3);
B = zeros(14,12);


B(1:3,1:3) = diag(ones(1,3)*timeStep);
B(4:6,4:6) = diag(ones(1,3)*timeStep);
B(7:9,1:3) = diag(ones(1,3)* timeStep^2/2);
B(7:9,7:9) = diag(ones(1,3)*timeStep);
B(10:12,4:6) = diag(ones(1,3)* timeStep^2/2);
B(10:12,10:12) = diag(ones(1,3)*timeStep);


newxk0_dot = subs(xk0_dot,[u v w p q r x y z phi theta ps]...
    ,[newu newv neww newp newq newr newx newy newz newphi newtheta newps]);
newxk0_dot = double(vpa(newxk0_dot));

Q = diag(std(B*(newxk0_dot*newxk0_dot')*B'))./1e6;

newCov = PHI*oldCov*PHI' + Q;


end

