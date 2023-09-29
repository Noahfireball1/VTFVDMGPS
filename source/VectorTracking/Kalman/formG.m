function G = formG(timeStep)

G = zeros(14,8);

G(1:6,1:6) = eye(6);
G(13:14,7:8) = eye(2);

G = G*timeStep;
end

