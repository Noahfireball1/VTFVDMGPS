function L = calcL(H,p_m,R)

S = H*p_m*H' + R;
L = p_m*H'*(S)^-1;

end

