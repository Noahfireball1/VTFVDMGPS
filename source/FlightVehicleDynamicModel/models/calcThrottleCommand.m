function throttle = calcThrottleCommand(refSpd,state)
%% Throttle Stick Command
aircraftSpeed = norm(state(1:3));

P = 0.75;
b = 0.574546090887085;

pGain = P*(b*refSpd - aircraftSpeed);

throttle = pGain;

if throttle > 1
    throttle = 1;
elseif throttle < 0
    throttle = 0;

end

