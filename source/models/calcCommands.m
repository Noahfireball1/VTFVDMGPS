function [bank,altitude,velocity,propEff,fuelEff,flapAngle,leftBrake,rightBrake,ailTrim,elevTrim,rudTrim] = calcCommands(target,position,currentVelocity)

ABC = [position(1:2)';target(:,1:2)];

[radius,~] = fit_circle_through_3_points(ABC);

velo = norm(currentVelocity);

bank = atand(velo^2/(9.81*radius));
altitude = -1*target(1,3);

velocity = 65;
propEff = 0.925;
fuelEff = 1;
flapAngle = 0;
leftBrake = 0;
rightBrake = 0;
ailTrim = 0;
elevTrim = 0;
rudTrim = 0;

end

