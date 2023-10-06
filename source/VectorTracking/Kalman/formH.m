function H = formH(unitVectors,LLA)
lat = LLA(1);
long = LLA(2);
alt = LLA(3);
count = 1;
e = 0.0818191910428; % eccentricity
a = 6378137.0; % equatorial radius [meters]

C_e_n = [-cos(long)*sin(lat) -sin(long)*sin(lat) cos(lat);...
    -sin(long)                     cos(long)           0      ;...
    -cos(long)*cos(lat)  -cos(lat)*sin(long) -sin(lat)];

R_N = (a*(1-e^2))/(1-e^2*sin(lat)^2)^(3/2); % [checked]
R_E = (a)/(1 - e^2*sin(lat)^2)^(1/2); % [checked]



for i = 1:2:2*length(unitVectors(1,:))

    uv_n = C_e_n*unitVectors(:,count); % unit vectors in the navigation frame
    uv_n_lla = [(R_N + alt)*uv_n(1);...
        (R_E + alt)*cos(lat)*uv_n(2);...
        -uv_n(3)];                        % unit vectors in the navigation frame with respect to LLA

    H(i,1:14) = [uv_n(1) uv_n(2) uv_n(3) zeros(1,3) zeros(1,3) zeros(1,3) 0 -1];
    H(i+1,1:14) = [zeros(1,3) zeros(1,3) uv_n_lla(1) uv_n_lla(2) uv_n_lla(3) zeros(1,3) -1 0];

    count = count + 1;
end

end

