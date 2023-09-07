classdef PropagateStates < handle
    %PROPAGATESTATES Summary of this class goes here
    %   Detailed explanation goes here

    methods
        function [obj,X_Dot] = PropagateStates(f_ib_b,m_ib_b,MOI,cg,invMassMatrix,mass,state,timeStep)

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

        end

        function covariance = calcCovariance(obj,MOI,cg,invMassMatrix,mass,timeStep,state)

            % J(1,1) = (1067922029527797*p)/45035996273704960 + (150584015030824269*r)/36028797018963968000 + (298977403082689*v)/1267650600228229401496703205376;
            % J(1,2) = 0;
            % J(1,3) = 0;
            % J(1,4) = -(479633360314957851*q)/562949953421312000 - (298977403082689*u)/1267650600228229401496703205376 - (7205759403792793197*w)/7205759403792793600;
            % J(1,5) = (150584015030824269*p)/36028797018963968000 - (985964771368110627*r)/1125899906842624000 + (7205759403792793197*v)/7205759403792793600;
            % J(1,6) = 0;
            % J(1,7) = 0;
            % J(1,8) = 0;
            % J(1,9) = 0;
            % J(1,10) = -(298977403082689*q)/1267650600228229401496703205376;
            % J(1,11) = (298977403082689*p)/1267650600228229401496703205376 + (7205759403792793197*r)/7205759403792793600;
            % J(1,12) = -(7205759403792793197*q)/7205759403792793600;
            %
            % J(2,1) = (82319691341016791801*q)/288230376151711744000 + (507060240091291686981746483975017*w)/507060240091291760598681282150400;
            % J(2,2) = 0;
            % J(2,3) = 0;
            % J(2,4) = (82319691341016791801*p)/288230376151711744000 - (30774959053168085347*r)/1475739525896764129280;
            % J(2,5) = -(30774959053168085347*q)/1475739525896764129280 - (507060240091291686981746483975017*u)/507060240091291760598681282150400;
            % J(2,6) = 0;
            % J(2,7) = 0;
            % J(2,8) = 0;
            % J(2,9) = 0;
            % J(2,10) = ;
            % J(2,11) = ;
            % J(2,12) = ;
            %
            % J(3,1) = ;
            % J(3,2) = ;
            % J(3,3) = ;
            % J(3,4) = ;
            % J(3,5) = ;
            % J(3,6) = ;
            % J(3,7) = ;
            % J(3,8) = ;
            % J(3,9) = ;
            % J(3,10) = ;
            % J(3,11) = ;
            % J(3,12) = ;
            %
            % J(4,1) = ;
            % J(4,2) = ;
            % J(4,3) = ;
            % J(4,4) = ;
            % J(4,5) = ;
            % J(4,6) = ;
            % J(4,7) = ;
            % J(4,8) = ;
            % J(4,9) = ;
            % J(4,10) = ;
            % J(4,11) = ;
            % J(4,12) = ;
            %
            % J(5,1) = ;
            % J(5,2) = ;
            % J(5,3) = ;
            % J(5,4) = ;
            % J(5,5) = ;
            % J(5,6) = ;
            % J(5,7) = ;
            % J(5,8) = ;
            % J(5,9) = ;
            % J(5,10) = ;
            % J(5,11) = ;
            % J(5,12) = ;
            %
            % J(6,1) = ;
            % J(6,2) = ;
            % J(6,3) = ;
            % J(6,4) = ;
            % J(6,5) = ;
            % J(6,6) = ;
            % J(6,7) = ;
            % J(6,8) = ;
            % J(6,9) = ;
            % J(6,10) = ;
            % J(6,11) = ;
            % J(6,12) = ;
            %
            % J(7,1) = ;
            % J(7,2) = ;
            % J(7,3) = ;
            % J(7,4) = ;
            % J(7,5) = ;
            % J(7,6) = ;
            % J(7,7) = ;
            % J(7,8) = ;
            % J(7,9) = ;
            % J(7,10) = ;
            % J(7,11) = ;
            % J(7,12) = ;
            %
            % J(8,1) = ;
            % J(8,2) = ;
            % J(8,3) = ;
            % J(8,4) = ;
            % J(8,5) = ;
            % J(8,6) = ;
            % J(8,7) = ;
            % J(8,8) = ;
            % J(8,9) = ;
            % J(8,10) = ;
            % J(8,11) = ;
            % J(8,12) = ;
            %
            % J(9,1) = ;
            % J(9,2) = ;
            % J(9,3) = ;
            % J(9,4) = ;
            % J(9,5) = ;
            % J(9,6) = ;
            % J(9,7) = ;
            % J(9,8) = ;
            % J(9,9) = ;
            % J(9,10) = ;
            % J(9,11) = ;
            % J(9,12) = ;
            %
            % J(10,1) = ;
            % J(10,2) = ;
            % J(10,3) = ;
            % J(10,4) = ;
            % J(10,5) = ;
            % J(10,6) = ;
            % J(10,7) = ;
            % J(10,8) = ;
            % J(10,9) = ;
            % J(10,10) = ;
            % J(10,11) = ;
            % J(10,12) = ;
            %
            % J(11,1) = ;
            % J(11,2) = ;
            % J(11,3) = ;
            % J(11,4) = ;
            % J(11,5) = ;
            % J(11,6) = ;
            % J(11,7) = ;
            % J(11,8) = ;
            % J(11,9) = ;
            % J(11,10) = ;
            % J(11,11) = ;
            % J(11,12) = ;
            %
            % J(12,1) = ;
            % J(12,2) = ;
            % J(12,3) = ;
            % J(12,4) = ;
            % J(12,5) = ;
            % J(12,6) = ;
            % J(12,7) = ;
            % J(12,8) = ;
            % J(12,9) = ;
            % J(12,10) = ;
            % J(12,11) = ;
            % J(12,12) = ;


            covariance = ones(12);


        end

    end
end

