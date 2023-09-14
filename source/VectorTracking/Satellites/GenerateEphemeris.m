classdef GenerateEphemeris < handle
    %GENERATEEPHEMERIS Summary of this class goes here
    %   Detailed explanation goes here
    properties (Access = public)
        rinex

    end

    properties (Access = private)
        ephYear;
        ephMonth;
        ephDay;
        login = [];
        netRCFilePath = [];

    end

    properties (Constant, Access = private)
        c = 299792458;
        wavelength = 299792458/1575.42e6;
        transmitFreq = 1575.42e6;

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
                    obj.rinex = rinexread(rinexFilePath);

                case ".Z"

            end

        end

        function rinexFile = loadRinexFile(obj,dir)

            yearSuffix = num2str(obj.ephYear);

            baseURL = "https://cddis.nasa.gov/archive/gnss/data/daily/";
            webURL = sprintf('%s/%s/%sn/',num2str(obj.ephYear),string(obj.julianDay),yearSuffix(3:end));
            fileURL = sprintf('AMC400USA_R_%s%s0000_01D_GN.rnx.gz',num2str(obj.ephYear),string(obj.julianDay));
            url = append(baseURL, webURL,fileURL);
            outputFilePath = fullfile(dir.rinex,fileURL);
            [path,unzippedFile,~] = fileparts(outputFilePath);

            if ~exist(append(path,filesep,unzippedFile),'file')

                % utilities.printText.options(1)

                obj.downloadRinex(url,dir,fileURL)

                obj.uncompressRinex(dir,fileURL)

                % utilities.printText.options(2)

            end

            rinexFile = append(path,filesep,unzippedFile);

        end

        function downloadRinex(obj,url,dir,fileURL)

            obj.checkPassword(dir);

            system(sprintf('curl -s -c %scookies.txt --ciphers DEFAULT@SECLEVEL=1 --netrc-file %s -L -O %s',dir.rinex,obj.netRCFilePath,url));

            system(sprintf('mv %s %s',fileURL,dir.rinex));

            system(sprintf('del %scookies.txt',dir.rinex));

        end

        function uncompressRinex(obj,dir,fileURL)

            switch obj.extension
                case ".rnx.gz"
                    gunzip(append(dir.rinex,fileURL));
                    system(sprintf('del %s%s',dir.rinex,fileURL));
                case ".gz"
                    uncompress(append(dir.rinex,fileURL));
            end

        end

        function checkPassword(obj,dir)

            if ~exist(append(dir.config,'CDDISLogin.txt'),"file")

                % utilities.printText.options(3)
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
    methods (Access = public)






    end
end

