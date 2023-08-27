function pedalPosn = calcPedalPosnCommand(atmos,state)

relativeVelocity = atmos.windVector - state(1:3);

sideSlip = asind(relativeVelocity(2)/norm(relativeVelocity));

pedalPosn = 0.25*(-sideSlip);

if pedalPosn > 1
    pedalPosn = 1;
elseif pedalPosn < 1
    pedalPosn = -1;
end

end