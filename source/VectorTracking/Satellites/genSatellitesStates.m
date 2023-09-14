function satStates = genSatellitesStates(time,date,dir)

% Pull Rinex for given date and time
rinex = GenerateEphemeris(date,dir);

timeArray = 0:1/50:time;

dayVector = repmat(day(date),length(timeArray),1);
monthVector = repmat(month(date),length(timeArray),1);
yearVector = repmat(year(date),length(timeArray),1);
hourVector = zeros(length(timeArray),1);
minuteVector = zeros(length(timeArray),1);
secondVector = timeArray';

dateVector = [yearVector monthVector dayVector hourVector minuteVector secondVector];

dateTimeVector = datetime(dateVector,'Format','d-MMM-y HH:mm:ss.SSS');
data = rinex.rinex.GPS;
[~,satIdx] = unique(data.SatelliteID);
data = data(satIdx,:);

for i = 1:length(timeArray)
    [satPos,satVel,satID] = gnssconstellation(dateTimeVector(i),data,"GNSSFileType","RINEX");

    satStates(:,:,i) = [satPos(:,1)';satVel(:,1)';satPos(:,2)';satVel(:,2)';satPos(:,3)';satVel(:,3)';zeros(1,length(satID))];
end
end

