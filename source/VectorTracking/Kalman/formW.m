function w = formW(Qd,S)
c = 299792458;
diagS = diag(S);
w = zeros(8,1);

w(1:3) = randn(3,1).*sqrt(diagS(1:3));
w(4:6) = randn(3,1).*sqrt(diagS(4:6));
w(7:8) = sqrt(Qd(13:14,13:14))*randn(2,1);

end

