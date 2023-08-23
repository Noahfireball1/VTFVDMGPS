classdef GenerateEphemeris < handle
    %GENERATEEPHEMERIS Summary of this class goes here
    %   Detailed explanation goes here
    properties (Access = public)
        eph;

    end

    properties (Access = private)
        ephYear = 1998;
        ephMonth = 12;
        ephDay = 22;
        login = [];
        netRCFilePath = [];

    end

    properties (Dependent, Access = private)
        julianDay;
        extension;
    end

    methods
        function obj = GenerateEphemeris(date,dir)

            obj.ephYear = year(date);
            obj.ephMonth = month(date);
            obj.ephDay = day(date);

            obj.loadEphemeris(dir)

        end

    end
    methods (Access = private)

        function loadEphemeris(obj,dir)

            % Load or Download Correct Rinex File
            rinexFilePath = obj.loadRinexFile(dir);

            % Parse Ephemeris from File
            switch obj.extension
                case ".rnx.gz"
                    ephTemp = rinexread(rinexFilePath);
                    ephTemp = ephTemp.GPS;

                case ".Z"

            end

            % Select only ephemeris at 0000 hrs
            timeIdxs = ephTemp.Time == datetime(obj.ephYear,obj.ephMonth,obj.ephDay,0,0,0);

            ephTemp = table2struct(ephTemp(timeIdxs,:));
            for i = 1:height(ephTemp)

                obj.eph.satelliteID(1,i) = ephTemp(i).SatelliteID;
                obj.eph.SVClockBias(1,i) = ephTemp(i).SVClockBias;
                obj.eph.SVClockDrift(1,i) = ephTemp(i).SVClockDrift;
                obj.eph.SVClockDriftRate(1,i) = ephTemp(i).SVClockDriftRate;
                obj.eph.IODE(1,i) = ephTemp(i).IODE;
                obj.eph.Crs(1,i) = ephTemp(i).Crs;
                obj.eph.Delta_n(1,i) = ephTemp(i).Delta_n;
                obj.eph.M0(1,i) = ephTemp(i).M0;
                obj.eph.Cuc(1,i) = ephTemp(i).Cuc;
                obj.eph.Eccentricity(1,i) = ephTemp(i).Eccentricity;
                obj.eph.Cus(1,i) = ephTemp(i).Cus;
                obj.eph.sqrtA(1,i) = ephTemp(i).sqrtA;
                obj.eph.Toe(1,i) = ephTemp(i).Toe;
                obj.eph.Cic(1,i) = ephTemp(i).Cic;
                obj.eph.OMEGA0(1,i) = ephTemp(i).OMEGA0;
                obj.eph.Cis(1,i) = ephTemp(i).Cis;
                obj.eph.i0(1,i) = ephTemp(i).i0;
                obj.eph.Crc(1,i) = ephTemp(i).Crc;
                obj.eph.omega(1,i) = ephTemp(i).omega;
                obj.eph.OMEGA_DOT(1,i) = ephTemp(i).OMEGA_DOT;
                obj.eph.IDOT(1,i) = ephTemp(i).IDOT;
                obj.eph.L2ChannelCodes(1,i) = ephTemp(i).L2ChannelCodes;
                obj.eph.GPSWeek(1,i) = ephTemp(i).GPSWeek;
                obj.eph.L2PDataFlag(1,i) = ephTemp(i).L2PDataFlag;
                obj.eph.SVAccuracy(1,i) = ephTemp(i).SVAccuracy;
                obj.eph.SVHealth(1,i) = ephTemp(i).SVHealth;
                obj.eph.TGD(1,i) = ephTemp(i).TGD;
                obj.eph.IODC(1,i) = ephTemp(i).IODC;
                obj.eph.TransmissionTime(1,i) = ephTemp(i).TransmissionTime;
                obj.eph.FitInterval(1,i) = ephTemp(i).FitInterval;
                obj.eph.BRDCOrbit7Spare3(1,i) = ephTemp(i).BRDCOrbit7Spare3;
                obj.eph.BRDCOrbit7Spare4(1,i) = ephTemp(i).BRDCOrbit7Spare4;


            end

        end

        function rinexFile = loadRinexFile(obj,dir)

            yearSuffix = num2str(obj.ephYear);

            baseURL = "https://cddis.nasa.gov/archive/gnss/data/daily/";
            webURL = sprintf('%s/%s/%sn/',num2str(obj.ephYear),string(obj.julianDay),yearSuffix(3:end));
            fileURL = sprintf('AMC400USA_R_%s%s0000_01D_GN.rnx.gz',num2str(obj.ephYear),string(obj.julianDay));
            url = append(baseURL, webURL,fileURL);
            outputFilePath = fullfile(dir.dataGPS,fileURL);
            [path,unzippedFile,~] = fileparts(outputFilePath);

            if ~exist(append(path,filesep,unzippedFile),'file')

                utilities.printText.options(1)

                obj.downloadRinex(url,dir,fileURL)

                obj.uncompressRinex(dir,fileURL)

                utilities.printText.options(2)

            end

            rinexFile = append(path,filesep,unzippedFile);

        end

        function downloadRinex(obj,url,dir,fileURL)

            obj.checkPassword(dir);

            system(sprintf('curl -s -c %scookies.txt --ciphers DEFAULT@SECLEVEL=1 --netrc-file %s -L -O %s',dir.dataGPS,obj.netRCFilePath,url));

            system(sprintf('mv %s %s',fileURL,dir.dataGPS));

            system(sprintf('del %scookies.txt',dir.dataGPS));

        end

        function uncompressRinex(obj,dir,fileURL)

            switch obj.extension
                case ".rnx.gz"
                    gunzip(append(dir.dataGPS,fileURL));
                    system(sprintf('del %s%s',dir.dataGPS,fileURL));
                case ".gz"
                    uncompress(append(dir.dataGPS,fileURL));
            end

        end

        function checkPassword(obj,dir)

            if ~exist(append(dir.config,'CDDISLogin.txt'),"file")

                utilities.printText.options(3)
                username = input('Username:','s');
                password = input('Password:','s');
                obj.createCDDISLoginFile(username,password,dir);
                obj.createNetRCFile(username,password,dir);
            elseif ~exist(append(dir.config,'.netrc'),"file")

                fid = fopen(append(dir.config,'CDDISLogin.txt'));

                loginFile = fscanf(fid,'%c');
                splitLogin = split(loginFile,"password:");

                username = splitLogin{1}(10:end-2);
                password = splitLogin{2}(1:end-2);

                obj.createNetRCFile(username,password,dir);
            end
        end

        function createCDDISLoginFile(~,username,password,dir)

            fid = fopen(append(dir.config,'CDDISLogin.txt'),"wt");
            fprintf(fid,'username:%s\npassword:%s\n',username,password);
            fclose(fid);

        end

        function createNetRCFile(obj,username,password,dir)

            fid = fopen(append(dir.config,'.netrc'),"wt");
            fprintf(fid,'machine urs.earthdata.nasa.gov login %s password %s',username,password);
            fclose(fid);

            obj.netRCFilePath = append(dir.config,'.netrc');

        end

    end

    methods

        function julianDay = get.julianDay(obj)
            julianDay = day(datetime(obj.ephYear,obj.ephMonth,obj.ephDay),"dayofyear");
        end

        function extension = get.extension(obj)
            if obj.ephYear > 2020
                extension = ".rnx.gz";
            else
                extension = ".gz";
            end
        end

    end
end

