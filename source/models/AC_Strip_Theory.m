function [Strip_Theory_Forces,Strip_Theory_Moments] = AC_Strip_Theory(rho,X,U,stripGeom)
STNumbers = load("DA40ST.mat");
STNumbers = STNumbers.STNumbers;

Vehicle = load("DA40.mat");
Vehicle = Vehicle.Vehicle;
%% Defining State Vector to Code Easier
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

V = norm([u;v;w]);
%% Defining Strip Theory Geometry from the excel and text files...
xp = stripGeom(:,1); % x location
yp = stripGeom(:,2); % y location
zp = stripGeom(:,3); % z location
S_refp = stripGeom(:,4); % Strip area
psi_s = stripGeom(:,5); % psi angle of each strip
theta_s = stripGeom(:,6); % theta angle of each strip
phi_s = stripGeom(:,7); % phi angle of each strip
cp = stripGeom(:,8); % quarter chord location of each strip
%% Defining Each strips velocity
u_p = u - r.*yp + q.*zp;
v_p = v + r.*xp - p.*zp;
w_p = w - q.*xp + p.*yp;



u_s = (cosd(theta_s).*cosd(psi_s)).*u_p + ...
      (cosd(theta_s).*sind(psi_s)).*v_p +...
       -sind(theta_s).*w_p;
v_s = (sind(phi_s).*sind(theta_s).*cosd(psi_s) - sind(psi_s).*cosd(phi_s)).*u_p +...
      (sind(phi_s).*sind(theta_s).*sind(psi_s) + cosd(psi_s).*cosd(phi_s)).*v_p +...
      (sind(phi_s).*cosd(theta_s)).*w_p; 
  
w_s = (cosd(phi_s).*sind(theta_s).*cosd(psi_s) + sind(phi_s).*sind(psi_s)).*u_p +...
      (cosd(phi_s).*sind(theta_s).*sind(psi_s) - sind(phi_s).*cosd(psi_s)).*v_p +...
      (cosd(phi_s).*cosd(theta_s)).*w_p;
  
  
%% Defining geometric angle of attack and dynamic pressure of each strip
AOA_p = atand(w_s./(u_s + 1e-4));
V_p = sqrt(u_s.^2 + 0*v_s.^2 + w_s.^2);
qbar_p = 1/2.*rho.*V_p.^2;
%% Control Surface Deflections
%% Aileron
Ail_Setting = U(2);
Ail_Rat = 0.3;

AOA_p(13:22) = AOA_p(13:22) - atand(Ail_Rat*sin(Ail_Setting)/(1 - Ail_Rat + Ail_Rat*cos(Ail_Setting)));
AOA_p(38:47) = AOA_p(38:47) + atand(Ail_Rat*sin(Ail_Setting)/(1 - Ail_Rat + Ail_Rat*cos(Ail_Setting)));
%% Elevator
Elev_Setting = U(1);
Elev_Rat = 0.3;

AOA_p(51:66) = AOA_p(51:66) + atand(Elev_Rat*sin(Elev_Setting)/(1 - Elev_Rat + Elev_Rat*cos(Elev_Setting)));
AOA_p(67:82) = AOA_p(67:82) + atand(Elev_Rat*sin(Elev_Setting)/(1 - Elev_Rat + Elev_Rat*cos(Elev_Setting)));
%% Rudder
Rud_Setting = U(3);
Rud_Rat = 0.3;

AOA_p(83:89) = AOA_p(83:89) + atand(Rud_Rat*sin(Rud_Setting)/(1 - Rud_Rat + Rud_Rat*cos(Rud_Setting)));
%% Flap
Flap_Setting = U(7);
Flap_Rat = 0.3;

AOA_p(4:12) = AOA_p(4:12) + atand(Flap_Rat*sind(Flap_Setting)/(1 - Flap_Rat + Flap_Rat*cosd(Flap_Setting)));
AOA_p(29:37) = AOA_p(29:37) + atand(Flap_Rat*sind(Flap_Setting)/(1 - Flap_Rat + Flap_Rat*cosd(Flap_Setting)));
%% Downwash script for calculating induced, then effective AOA
% 1-25: Left wing, root to tip
LWIndices = 1:25;
% 26-50: Right wing, root to tip
RWIndices = 26:50;
% 51-66: Left horizontal stab, root to tip
LHTIndices = 51:66;
% 67-82: Right horizontal stab, root to tip
RHTIndices = 67:82;
% 83-89: Vertical stab, root to tip
VTIndices = 83:89;

% NOTE AOA_p: the "point" or "geometric" angles of attack that you have
% computed for each strip

% Get induced angles of attack
AOAi_LeftWing = CalcInducedAOA(Vehicle.Aero.StripROM.LeftWing,AOA_p(LWIndices),V_p(LWIndices));
AOAi_RightWing = CalcInducedAOA(Vehicle.Aero.StripROM.RightWing,AOA_p(RWIndices),V_p(RWIndices));

AOAi_LeftHT = CalcInducedAOA(Vehicle.Aero.StripROM.LeftHT,AOA_p(LHTIndices),V_p(LHTIndices));
AOAi_RightHT = CalcInducedAOA(Vehicle.Aero.StripROM.RightHT,AOA_p(RHTIndices),V_p(RHTIndices));

AOAi_VT = CalcInducedAOA(Vehicle.Aero.StripROM.VT,AOA_p(VTIndices),V_p(VTIndices));

AOA_ind = [AOAi_LeftWing; AOAi_RightWing; AOAi_LeftHT; AOAi_RightHT; AOAi_VT];

% compute the effective angle of attack --> Effective AOA = Geom AOA - Induced AOA
AOA_eff = AOA_p - AOA_ind;

% ^^^^ These are the effective angles of attack seen by each strip after
% accounting for downwash. Use THESE to query the sectional aero
% characteristics of the strips. 

%% Importing Aerodynamic Coefficients
%% Right Wing
rWing_in_CL = STNumbers.clWing(AOA_eff(26:28));
rWing_in_CD = STNumbers.cdWing(AOA_eff(26:28));
rWing_in_Cm = STNumbers.cmyWing(AOA_eff(26:28));
rWing_Flap_CL = STNumbers.clWing(AOA_eff(29:37));
rWing_Flap_CD = STNumbers.cdWing(AOA_eff(29:37));
rWing_Flap_Cm = STNumbers.cmyWing(AOA_eff(29:37));
rWing_Aileron_CL = STNumbers.clWing(AOA_eff(38:47));
rWing_Aileron_CD = STNumbers.cdWing(AOA_eff(38:47));
rWing_Aileron_Cm = STNumbers.cmyWing(AOA_eff(38:47));
rWing_out_CL = STNumbers.clWing(AOA_eff(48:50));
rWing_out_CD = STNumbers.cdWing(AOA_eff(48:50));
rWing_out_Cm = STNumbers.cmyWing(AOA_eff(48:50));
%% Left Wing
lWing_in_CL = STNumbers.clWing(AOA_eff(1:3));
lWing_in_CD = STNumbers.cdWing(AOA_eff(1:3));
lWing_in_Cm = STNumbers.cmyWing(AOA_eff(1:3));
lWing_Flap_CL = STNumbers.clWing(AOA_eff(4:12));
lWing_Flap_CD = STNumbers.cdWing(AOA_eff(4:12));
lWing_Flap_Cm = STNumbers.cmyWing(AOA_eff(4:12));
lWing_Aileron_CL = STNumbers.clWing(AOA_eff(13:22));
lWing_Aileron_CD = STNumbers.cdWing(AOA_eff(13:22));
lWing_Aileron_Cm = STNumbers.cmyWing(AOA_eff(13:22));
lWing_out_CL = STNumbers.clWing(AOA_eff(23:25));
lWing_out_CD = STNumbers.cdWing(AOA_eff(23:25));
lWing_out_Cm = STNumbers.cmyWing(AOA_eff(23:25));
%% Right Horizantal Tail
rTail_in_CL = STNumbers.clHor(AOA_eff(67:82));
rTail_in_CD = STNumbers.cdHor(AOA_eff(67:82));
rTail_in_Cm = STNumbers.cmyHor(AOA_eff(67:82));
%% Left Horizantal Tail
lTail_in_CL = STNumbers.clHor(AOA_eff(51:66));
lTail_in_CD = STNumbers.cdHor(AOA_eff(51:66));
lTail_in_Cm = STNumbers.cmyHor(AOA_eff(51:66));
%% Vertical Tail
vTail_in_CL = STNumbers.clVert(AOA_eff(83:89));
vTail_in_CD = STNumbers.cdVert(AOA_eff(83:89));
vTail_in_Cm = STNumbers.cmyVert(AOA_eff(83:89));
%% Assimilating all of these Aerodynamci Coeffs into singular matrices...
%% Lift Coefficient
CL = [lWing_in_CL;lWing_Flap_CL;lWing_Aileron_CL;lWing_out_CL;rWing_in_CL;rWing_Flap_CL;rWing_Aileron_CL;rWing_out_CL;lTail_in_CL;rTail_in_CL;vTail_in_CL];
%% Drag Coefficient
CD = [lWing_in_CD;lWing_Flap_CD;lWing_Aileron_CD;lWing_out_CD;rWing_in_CD;rWing_Flap_CD;rWing_Aileron_CD;rWing_out_CD;lTail_in_CD;rTail_in_CD;vTail_in_CD];
%% Moment Coefficient
Cm = [lWing_in_Cm;lWing_Flap_Cm;lWing_Aileron_Cm;lWing_out_Cm;rWing_in_Cm;rWing_Flap_Cm;rWing_Aileron_Cm;rWing_out_Cm;lTail_in_Cm;rTail_in_Cm;vTail_in_Cm];
%% Lifting Force at Each Strip
Ls = qbar_p.*S_refp.*CL; 
%% Drag Force at Each Strip
Ds = qbar_p.*S_refp.*CD;
%% Moment at the c/4 of each strip
Ms_c4 = qbar_p.*S_refp.*cp.*Cm; 
%% Summing the lift and drag of each individual strip
Fx_s = Ls.*sind(AOA_eff) - Ds.*cosd(AOA_eff);
Fz_s = -Ls.*cosd(AOA_eff) - Ds.*sind(AOA_eff);
My_s = Ms_c4;
%% Doing a 321 Rotation to put forces into body axis
Fx_b = (cosd(theta_s).*cosd(psi_s)).*Fx_s + (sind(phi_s).*sind(theta_s).*cosd(psi_s) - sind(psi_s).*cosd(phi_s)).*0 + (cosd(phi_s).*sind(theta_s).*cosd(psi_s) + sind(phi_s).*sind(psi_s)).*Fz_s; 
Fy_b = (cosd(theta_s).*sind(psi_s)).*Fx_s + (sind(phi_s).*sind(theta_s).*sind(psi_s) + cosd(psi_s).*cosd(phi_s)).*0+ (cosd(phi_s).*sind(theta_s).*sind(psi_s) - sind(phi_s).*cosd(psi_s)).*Fz_s; 
Fz_b = -sind(theta_s).*Fx_s + (sind(phi_s).*cosd(theta_s)).*0 + (cosd(phi_s).*cosd(theta_s)).*Fz_s;
%% Doing the same for the moments...
Mx_c4s_b = (cosd(theta_s).*cosd(psi_s)).*0 + (sind(phi_s).*sind(theta_s).*cosd(psi_s) - sind(psi_s).*cosd(phi_s)).*My_s + (cosd(phi_s).*sind(theta_s).*cosd(psi_s) + sind(phi_s).*sind(psi_s)).*0;        
My_c4s_b = (cosd(theta_s).*sind(psi_s)).*0 + (sind(phi_s).*sind(theta_s).*sind(psi_s) + cosd(psi_s).*cosd(phi_s)).*My_s+ (cosd(phi_s).*sind(theta_s).*sind(psi_s) - sind(phi_s).*cosd(psi_s)).*0;         
Mz_c4s_b = -sind(theta_s).*0 + (sind(phi_s).*cosd(theta_s)).*My_s + (cosd(phi_s).*cosd(theta_s)).*0;
M_c4s_b = [sum(Mx_c4s_b);sum(My_c4s_b);sum(Mz_c4s_b)];

F_lifting = [sum(Fx_b); sum(Fy_b); sum(Fz_b)];
M_lifting = [sum(-zp.*Fy_b + yp.*Fz_b);sum(zp.*Fx_b - xp.*Fz_b);sum(-yp.*Fx_b + xp.*Fy_b)] + M_c4s_b;
%% Repeating the Same thing for non-lifting section of the DA-40...
c_ref = 1.067;
S_ref = 12.91;
q_bar = 1/2*rho*(V^2);
alpha = atand(w/(u + 1e-4));
beta = asind(v/(V + 1e-4));

Cx = STNumbers.cx(alpha,abs(beta));
Cy = STNumbers.cy(alpha,abs(beta))*sign(beta);
Cz = STNumbers.cz(alpha,abs(beta));
Cmx = STNumbers.cmx(alpha,abs(beta))*sign(beta);
Cmy = STNumbers.cmy(alpha,abs(beta));
Cmz = STNumbers.cmz(alpha,abs(beta))*sign(beta);

nonlift_Fx = q_bar*S_ref*Cx;
nonlift_Fy = q_bar*S_ref*Cy;
nonlift_Fz = q_bar*S_ref*Cz;
nonlift_Fmx = q_bar*S_ref*Cmx*c_ref;
nonlift_Fmy = q_bar*S_ref*Cmy*c_ref;
nonlift_Fmz = q_bar*S_ref*Cmz*c_ref;

F_nonlift = [nonlift_Fx;nonlift_Fy;nonlift_Fz];
M_nonlift = [nonlift_Fmx;nonlift_Fmy;nonlift_Fmz];
%% Calculating the Dampening characteristic of the control surfaces...
b = 11.32;
c_bar = 1.067;
%% Had to put this in here to stop throwing errors...
if V < 0.5
    p_hat = 0;
    q_hat = 0;
    r_hat = 0;
else
    p_hat = p*b/(2*V);
    q_hat = q*c_bar/(2*V);
    r_hat = r*b/(2*V);
end
%% Using interpolation functions for all the hat coeffs...
CLq_ref = [0.527 1.8001 1.0984 0.4558];
Cmq_ref = [0.195 0.1973 0.1989 0.1997];
Cyr_ref = [-0.005 0.0535 0.0007 0.0097];
Clr_ref = [-0.0004 -0.0003 -0.0002 -0.0001];
Cnr_ref = [-0.0013 -0.0013 -0.0012 -0.0011];
AOA_ref = [0,2,3,5];
Cyp_ref = [0.0269 0.0127 0.037 -0.0387];
Cyp = interp1(AOA_ref,Cyp_ref,alpha, 'linear', 'extrap');

Clp_ref = [0 0 0 0];
Clp = interp1(AOA_ref,Clp_ref,alpha, 'linear', 'extrap'); 

Cnp_ref = [0.0008 0.0007 0.0007 0.0007];
Cnp = interp1(AOA_ref,Cnp_ref,alpha, 'linear', 'extrap'); 
CLq = interp1(AOA_ref,CLq_ref,alpha, 'linear', 'extrap'); 
Cmq = interp1(AOA_ref,Cmq_ref,alpha, 'linear', 'extrap');
Cyr = interp1(AOA_ref,Cyr_ref,alpha, 'linear', 'extrap'); 
Clr = interp1(AOA_ref,Clr_ref,alpha, 'linear', 'extrap');
Cnr = interp1(AOA_ref,Cnr_ref,alpha, 'linear', 'extrap'); 

CD_damp = 0;
CL_damp = CLq*q_hat;
CY_damp = Cyp*p_hat + Cyr*r_hat;

Cm_damp = Cmq*q_hat;
Cl_damp = Clp*p_hat + Clr*r_hat;
Cn_damp = Cnp*p_hat + Cnr*r_hat;

F_Damping = [CD_damp;CY_damp;CL_damp].*q_bar.*S_ref;
M_Damping = [Cl_damp*b;Cm_damp*c_bar;Cn_damp*b].*q_bar.*S_ref;
%% Adding Each Component of the sections and forces and moments together...
Strip_Theory_Forces = F_lifting + F_nonlift + F_Damping;
Strip_Theory_Moments = M_lifting + M_nonlift + M_Damping;

stopper = 1;
end