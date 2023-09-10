
function AOAi = CalcInducedAOA(StripROM,AOA,StripVel)

nStrips = height(StripROM(:,1));

Strips.Y      = StripROM(:,1);
Strips.Length = StripROM(:,2);
Strips.c0     = StripROM(:,3);
Strips.dCdA   = StripROM(:,4);

% Calculating the dA wrt to the reference -5 degree
dA = AOA + 5;                           % 5 AoA (specified angle) - (-5 AoA (baseline) )

%% COMPUTE CIRCULATION GRADIENTS FOR EACH STRIP "dCdY"

% dCdY = zeros(nStrips,1);

c1 = (Strips.c0(1:nStrips - 1) + Strips.dCdA(1:nStrips - 1).*dA(1:nStrips - 1)).*StripVel(1:nStrips - 1);
c2 = (Strips.c0(2:nStrips) + Strips.dCdA(2:nStrips).*dA(2:nStrips)).*StripVel(2:nStrips);
dy = Strips.Length(1:nStrips - 1);
dCdY = (c2 - c1)./dy;
dCdY(nStrips) = dCdY(nStrips - 1);


%% COMPUTE DOWNWASH VELOCITY FOR EACH STRIP

% Downwash = zeros(nStrips,1);
y0 = Strips.Y;

x = ((1./(Strips.Y' - Strips.Y)).* (dCdY.* Strips.Length))./(4*pi);
x(isinf(x)|isnan(x)) = 0;
Downwash = sum(x)';

% % For Root and Tip Correction
%

%     Compute the circulation values of the start and end strips.
% y1 = Strips.Y(1);
y2 = Strips.Y(nStrips);

%     compute the circulation values of the start and end strips.
% c0 = Strips.c0(1);
% dCdA = Strips.dCdA(1);

% c1 = (c0 + dCdA*dA).*StripVel;

c0 = Strips.c0(nStrips);
dCdA = Strips.dCdA(nStrips);

c2 = (c0 + dCdA*dA(nStrips)).*StripVel(nStrips);

% % Starting strip vortex contribution
%     (Note: don't add this if the
%     starting strip is connected to some other non-lifting component like fuselage, pod nacelle etc.)
%
% Downwash(1) = Downwash(1) + ((2./(Strips.Length(1)))*c1)/(4*pi);
%
% Downwash(2:end) = Downwash(2:end) + ((1./((y0(2:end) - y1)))*c1)/(4*pi);




% % Ending strip vortex contribution
%     (Note: don't add this if the end strip
%      is connected to some other non-lifting component like fuselage, pod nacelle etc.)

% This tip correction is heavily hardcoded for this model.
if nStrips > 7      % No tip correction for Vertical Tail
    Downwash(1:nStrips-1) = Downwash(1:nStrips-1) - ((1./((y0(1:nStrips-1) - y2)))*c2)/(4*pi);
    
    if sign(StripROM(nStrips,1)) == -1
        Downwash(nStrips) = Downwash(nStrips) - ((2./(Strips.Length(nStrips)))*c2)/(4*pi);
    else
        Downwash(nStrips) = Downwash(nStrips) + ((2./(Strips.Length(nStrips)))*c2)/(4*pi);
    end
end

% Total Downwash After adding the root/tip correction
Strips.Downwash = Downwash;



%% COMPUTE INDUCED LIFT & DRAG FROM ENTIRE COMPONENT

% circulation = (Strips.c0 + Strips.dCdA*dA).*StripVel;
dw = Strips.Downwash;
% dy = Strips.Length;

% Lift = rho.*StripVel.*circulation.*dy;                  % Lift per Strip
% InducedDrag = rho.*dw.*circulation.*dy;                % Induced Drag per Strip


% % Induced AOA

AOAi = atand((dw)./StripVel);


end