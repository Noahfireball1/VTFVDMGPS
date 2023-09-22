function fllError = calcFLLError(lastIP, IP, lastQP, QP)
pdiTime = 1/50;

fllError = (IP.*lastQP - lastIP.*QP)/(pi*pdiTime);
end