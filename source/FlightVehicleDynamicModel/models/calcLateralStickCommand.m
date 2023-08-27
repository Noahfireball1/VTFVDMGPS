function latStick = calcLateralStickCommand(desiredCourse,state)
R2D = 180/pi;
%% Lateral Stick Command

referenceBank = (desiredCourse - state(12))*R2D;

if referenceBank > 30
    referenceBank = 30;
elseif referenceBank < -30
    referenceBank = -30;
end

P = 0.107891692261891;
b = 0.548332573243903;

pGain = P*(b*referenceBank - state(10)*R2D);


latStick = pGain;

if latStick > 1
    latStick = 1;
elseif latStick < -1
    latStick = -1;

end

