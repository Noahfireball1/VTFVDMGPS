function fllError = calcFLLError(lastIP, IP, lastQP, QP)
pdiTime = 1/50;

cross = lastIP.*QP - IP.*lastQP;
dot = lastIP.*IP + lastQP.*QP;
fllError = atan2(cross,dot)./(pdiTime*2*pi);
end