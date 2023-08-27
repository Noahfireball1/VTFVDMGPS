function longStick = calcLongitudinalStickCommand(selWaypoint,state)
R2D = 180/pi;
%% Pitch Hold Command
refAlt = -selWaypoint(3);

P = 0.392178102036374;
b = 1;

pGain = P*(b*refAlt + state(9));

referencePitch = pGain;

if referencePitch > 15
    referencePitch = 15;
elseif referencePitch < -15
    referencePitch = -15;
end

%% Longitudinal Stick Command

P = 0.4;
b = 0.548332573243903;

pGain = P*(b*referencePitch - state(11)*R2D);

longStick = pGain;

if longStick > 1
    longStick = 1;
elseif longStick < -1
    longStick = -1;

end

