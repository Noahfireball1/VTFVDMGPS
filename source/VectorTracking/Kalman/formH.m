function H = formH(unitVectors,LLA)
lat = LLA(1);
long = LLA(2);
alt = LLA(3);

e = 0.0818191910428; % eccentricity
R_0 = 6378137.0; % equatorial radius [meters]

% Matrix Rotating from ECEF to Navigation Frame (Groves Eq. 2.150)
C_e_n = [-cos(long)*sin(lat) -sin(long)*sin(lat) cos(lat);...
    -sin(long)                     cos(long)           0      ;...
    -cos(long)*cos(lat)  -cos(lat)*sin(long) -sin(lat)];

% Meridian Radius of Curvature (Groves Eq. 2.105)
R_N = (R_0*(1-e^2))/(1 - (e^2)*sin(lat)^2)^(3/2);

% Tranverse Radius of Curvature (Groves Eq. 2.106)
R_E = (R_0)/(sqrt(1 - (e^2)*sin(lat)^2));


count = 1;
for i = 1:2:2*length(unitVectors(:,1))

    % Rotating the Unit Vectors for the Carrier Residual Measurements (Groves Eq. 8.39)
    uv_n = C_e_n*unitVectors(count,:)';

    % Rotating the Unit Vectors for the Pseudorange Residual Measurements (Groves Eq. 14.128)
    uv_n_lla = [(R_N + alt)*uv_n(1);...
        (R_E + alt)*cos(lat)*uv_n(2);...
        -uv_n(3)];     

    H(i,1:14) = [uv_n(1) uv_n(2) uv_n(3) zeros(1,3) zeros(1,3) zeros(1,3) 0 1];
    H(i+1,1:14) = [zeros(1,3) zeros(1,3) uv_n_lla(1) uv_n_lla(2) uv_n_lla(3) zeros(1,3) 1 0];

    count = count + 1;
end

end

