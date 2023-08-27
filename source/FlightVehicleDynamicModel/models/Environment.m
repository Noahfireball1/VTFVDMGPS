classdef Environment < handle
    %ENVIRONMENT Summary of this class goes here
    %   Detailed explanation goes here

    properties
        temperature
        pressure
        density
        vSound
        windVector
    end

    methods
        function obj = Environment(type,lat,long,alt,varAlt,varWind)

            switch upper(type)
                case 'TRUTH'
                    obj.calcAtmos(alt,[])
                    obj.calcWind(lat,long,alt,[])

                case 'NOISE'
                    obj.calcAtmos(alt,varAlt)
                    obj.calcWind([],[],[],varWind)

            end

        end

        function calcAtmos(obj,alt,varAlt)
            if isempty(varAlt)
                [T,a,P,rho] = atmoscoesa(alt);

                obj.temperature = T;
                obj.vSound = a;
                obj.pressure = P;
                obj.density = rho;

            else
                [T,a,P,rho] = atmoscoesa(alt + randn(1)*varAlt);

                obj.temperature = T;
                obj.vSound = a;
                obj.pressure = P;
                obj.density = rho;

            end

        end

        function calcWind(obj,lat,long,alt,varWind)
            if isempty(varWind)
                xyWind = atmoshwm(lat,long,alt);
                obj.windVector = [xyWind 0];

            else
                obj.windVector = randn([1,3])*varWind;

            end


        end
    end
end

