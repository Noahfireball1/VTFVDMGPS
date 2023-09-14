function H = formH(unitVectors)
count = 1;

for i = 1:2:2*length(unitVectors(1,:))
    ux = unitVectors(1,count);
    uy = unitVectors(2,count);
    uz = unitVectors(3,count);

    H(i,1:14) = [ux uy uz 0 0 0 0 0 0 0 0 0 -1 0];
    H(i+1,1:14) = [0 0 0 0 0 0 ux uy uz 0 0 0 0 -1];

    count = count + 1;
end
end

