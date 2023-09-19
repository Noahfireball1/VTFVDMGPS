function satStates = genSatellitesStates(time,year,month,day,rinexFilePath)

hourVector = 0;
minuteVector = 0;
secondVector = time;

dateVector = [year month day hourVector minuteVector secondVector];

currentDateTime = datetime(dateVector,'Format','d-MMM-y HH:mm:ss.SSS');
rinex = rinexread(rinexFilePath);
rinex = rinex.GPS;
[~,satIdx] = unique(rinex.SatelliteID);
data = rinex(satIdx,:);

[satPos,satVel,satID] = gnssconstellation(currentDateTime,data,"GNSSFileType","RINEX");

satStates = [satPos(:,1)';satVel(:,1)';satPos(:,2)';satVel(:,2)';satPos(:,3)';satVel(:,3)';zeros(1,length(satID))];
end


