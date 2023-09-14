function R = formR(varPsr,varCarr)
count = 1;
Rtmp = [];

for i = 1:2:2*length(varPsr)

    Rtmp = [Rtmp varCarr(count) varPsr(count)];

    count = count + 1;
end

R = diag(Rtmp);
end

