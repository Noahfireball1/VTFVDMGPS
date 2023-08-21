classdef VectorTracking < handle
    
    properties
        date;
        ephemeris;
        trueTrajectory;
        FVDMImu;
        
    end
    
    methods
        function obj = VectorTracking(configFilePath,dirs)

            % Read Selected Configuration File
            config = utilities.yamlmatlab.ReadYaml(configFilePath);

            % Parse Yaml File
            obj.parseYaml(config);

            % Download Ephemeris File from Internet
            obj.ephemeris = GenerateEphemeris(obj.date,dirs).eph;
                


        end
        
        function parseYaml(obj,config)

            gen = config.general;

            obj.date = datetime(gen.year,gen.month,gen.day);

            FVDMvariables = load(append('data',filesep,'FVDM',filesep,config.FVDMName));
            
            




        end
    end
end

