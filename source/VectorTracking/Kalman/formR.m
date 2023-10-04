function R = formR(varPsr,varCarr,varAng)
count = 1;
Rtmp = [];

for i = 1:2:2*length(varPsr)

    Rtmp = [Rtmp varCarr(count) varPsr(count)];

    count = count + 1;
end

R = diag([Rtmp varAng varAng(1) varAng(4) varAng(4) varAng(4)]);
end

