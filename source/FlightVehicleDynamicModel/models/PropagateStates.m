classdef PropagateStates < handle
    %PROPAGATESTATES Summary of this class goes here
    %   Detailed explanation goes here

    methods
        function [obj,state,covariance] = PropagateStates(f_ib_b,m_ib_b,MOI,cg,invMassMatrix,mass,state,timeStep)

            u = state(1);
            v = state(2);
            w = state(3);
            p = state(4);
            q = state(5);
            r = state(6);
            phi = state(10);
            theta = state(11);
            psi = state(12);

            Xcg = cg(1);
            Ycg = cg(2);
            Zcg = cg(3);

            CGssm = [0 -Zcg Ycg;Zcg 0 -Xcg;-Ycg Xcg 0];

            lambdaX = f_ib_b(1)/mass + r*v - q*w - 0 - 2*q*0 + 2*r*0 + Xcg*(q^2 + r^2) - Ycg*p*q - Zcg*p*r;
            lambdaY = f_ib_b(2)/mass + p*w - r*u - 0 - 2*r*0 + 2*p*0 - Xcg*p*q + Ycg*(p^2 +r^2) - Zcg*q*r;
            lambdaZ = f_ib_b(3)/mass + q*u - p*v - 0 - 2*p*0 + 2*q*0 - Xcg*p*r - Ycg*q*r + Zcg*(p^2 + q^2);

            MomentVector = [p;q;r];

            MomentVectorssm = [0 -r q;r 0 -p;-q p 0];

            muVector = [m_ib_b(1);m_ib_b(2);m_ib_b(3)] - (MomentVectorssm*MOI*MomentVector) - mass*CGssm*MomentVectorssm*[u;v;w];
            LAMEparams = [lambdaX;lambdaY;lambdaZ;muVector];
            X_Dot(1:6) = invMassMatrix*LAMEparams;

            xdot = cos(theta)*cos(psi)*u + sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi)*v + cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi)*w;
            ydot = cos(theta)*sin(psi)*u + sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi)*v + cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi)*w;
            zdot = -sin(theta)*u + sin(phi)*cos(theta)*v + cos(phi)*cos(theta)*w;
            phidot = p + q*sin(phi)*tan(theta) + r*cos(phi)*tan(theta);
            thdot = q*cos(phi) - r*sin(phi);
            psidot = q*sin(phi)/cos(theta) + r*cos(phi)/cos(theta);
            X_Dot(7:12) = [xdot;ydot;zdot;phidot;thdot;psidot];

            state = X_Dot'*timeStep + state;
        end

    end
end

