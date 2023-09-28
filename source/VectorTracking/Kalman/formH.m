function H = formH(unitVectors,LLA)
count = 1;
e = 0.0818191910428; % eccentricity
a = 6378137.0; % equatorial radius [meters]

 C_e_n = [-cos(LLA(2))*sin(LLA(1)) -sin(LLA(2))*sin(LLA(1)) cos(LLA(1));...
          sin(LLA(2))                     cos(LLA(2))           0      ;...
          cos(LLA(2))*cos(LLA(1))  cos(LLA(1))*sin(LLA(2)) sin(LLA(1))];

R_N = (a*(1-e^2))/(1-e^2*sin(LLA(2))^2)^(3/2); % [checked]
R_E = (a)/(1 - e^2*sin(LLA(2))^2)^(1/2); % [checked]



for i = 1:2:2*length(unitVectors(1,:))

    uv_n = C_e_n*unitVectors(:,count); % unit vectors in the navigation frame
    uv_n_lla = [(R_N + LLA(3))*uv_n(1);(R_E + LLA(3))*cos(LLA(2))*uv_n(2);-uv_n(3)];                        % unit vectors in the navigation frame with respect to LLA

    H(i,1:14) = [uv_n(1) uv_n(2) uv_n(3) 0 0 0 0 0 0 0 0 0 0 1];
    H(i+1,1:14) = [0 0 0 0 0 0 uv_n_lla(1) uv_n_lla(2) uv_n_lla(3) 0 0 0 1 0];

    count = count + 1;
end

end

