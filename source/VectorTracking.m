classdef VectorTracking < handle

    properties
        date;
        timeSteps;
        ephemeris;
        FVDM;
        currentSeed;
        startTime = 0;
        endTime;
        timeStep;


    end

    properties (Access = private)

        

    end

    methods
        function obj = VectorTracking(configFilePath,dirs)

            % Read Selected Configuration File
            config = utilities.yamlmatlab.ReadYaml(configFilePath);

            % Parse Yaml File
            obj.parseYaml(config,dirs);

            % Download Ephemeris File from Internet
            obj.ephemeris = GenerateEphemeris(obj.date,dirs).eph;

        end

        function process(obj)

            obj.currentSeed = randi(1e6);

            % Start Simulation
            results = sim('source\FlightVehicleDynamicModel\DA40_Flight_Model.slx');

        end


    end

    methods (Access = private)
        function parseYaml(obj,config,dirs)

            gen = config.general;
            obj.endTime = config.aircraft.time;
            obj.timeStep = 1/config.aircraft.frequency;

            obj.date = datetime(gen.year,gen.month,gen.day);
            obj.FVDM.dayOfYear = day(obj.date,"dayofyear");

            obj.FVDM.wind.velocity = config.noise.wind.velocity;
            obj.FVDM.wind.velocityTS = config.noise.wind.velocityTS;
            obj.FVDM.wind.direction = config.noise.wind.direction;
            obj.FVDM.wind.directionTS = config.noise.wind.directionTS;

            obj.FVDM.gravity.force = str2num(cell2mat(config.noise.gravity.force));
            obj.FVDM.gravity.forceTS = config.noise.gravity.forceTS;
            obj.FVDM.gravity.moment = str2num(cell2mat(config.noise.gravity.moment));
            obj.FVDM.gravity.momentTS = config.noise.gravity.momentTS;

            obj.FVDM.aero.force = str2num(cell2mat(config.noise.aero.force));
            obj.FVDM.aero.forceTS = config.noise.aero.forceTS;
            obj.FVDM.aero.moment = str2num(cell2mat(config.noise.aero.moment));
            obj.FVDM.aero.momentTS = config.noise.aero.momentTS;

            obj.FVDM.engine.force = str2num(cell2mat(config.noise.engine.force));
            obj.FVDM.engine.forceTS = config.noise.engine.forceTS;
            obj.FVDM.engine.moment = str2num(cell2mat(config.noise.engine.moment));
            obj.FVDM.engine.momentTS = config.noise.engine.momentTS;

            obj.FVDM.initialStates = [config.aircraft.initialState.u;...
                config.aircraft.initialState.v;...
                config.aircraft.initialState.w;...
                config.aircraft.initialState.p;...
                config.aircraft.initialState.q;...
                config.aircraft.initialState.r;...
                config.aircraft.initialState.x;...
                config.aircraft.initialState.y;...
                config.aircraft.initialState.z;...
                config.aircraft.initialState.phi;...
                str2num(config.aircraft.initialState.theta);...
                str2num(config.aircraft.initialState.psi)];

            Vehicle = load("DA40.mat");
            BSFC_LUT = load("DA40ENGINE.mat");
            PowerFactor_LUT = load("DA40PF.mat");
            STNumbers = load("DA40ST.mat");
            STGeometry = load("DA40STGEOM.mat");
            selWaypoints = load(sprintf('%s.mat',config.aircraft.waypoints));
            




            obj.FVDM.Vehicle = Vehicle.Vehicle;
            obj.FVDM.BSFC_LUT = BSFC_LUT.BSFC_LUT;
            obj.FVDM.PowerFactor_LUT = PowerFactor_LUT.PowerFactor_LUT;
            obj.FVDM.STNumbers = STNumbers.STNumbers;
            obj.FVDM.STGeometry = STGeometry.ST_Geometry;
            obj.FVDM.waypoints = selWaypoints.waypoints;
            obj.FVDM.lookaheadDist = config.aircraft.lookaheadDistance;
            obj.FVDM.refLL = selWaypoints.refLL;




        end

    end
end

