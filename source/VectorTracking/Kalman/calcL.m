function L = calcL(H,p_m,R)

S = H*p_m*H' + R;
L = p_m*H'*(S)^-1;

if isnan(L)
    L = zeros(size(L));
end

end

