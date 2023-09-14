function p_m = updateCovariance(p_m,L,H)

p_m = (eye(size(p_m)) - L*H)*p_m;

end

