function varAng = calcVarAng(CN0)
avgCN0 = mean(CN0);
varAng = [1/avgCN0.^2 1/avgCN0.^2 1/avgCN0.^2 1/abs(sqrt(avgCN0.^2)) 1/abs(sqrt(avgCN0.^2)) 1/abs(sqrt(avgCN0.^2))];
end

