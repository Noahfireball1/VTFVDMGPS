function Qd = formQ(w,Gamma)

Qd = diag(mean((Gamma*w*w'*Gamma')));
end

