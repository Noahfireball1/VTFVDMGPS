classdef VectorTrackingPlotting < handle
    %VECTORTRACKINGPLOTTING Summary of this class goes here
    %   Detailed explanation goes here

    properties
        outputs
        numMonteCarlos
    end

    methods
        function obj = VectorTrackingPlotting(outputFilePath)

            loadedFile = load(outputFilePath);
            obj.outputs = loadedFile.run;

            obj.numMonteCarlos = 1;

            obj.plotStates();

            obj.plotErrors();

            obj.plotCovariance();

            obj.plotResiduals();




        end

        function plotStates(obj)
            stateTimeArray = obj.outputs.tout(1):1/50:obj.outputs.tout(end);
            figure
            tiledlayout(3,1)
            nexttile
            hold on
            ylabel('North Velocity [m/s]')
            plot(stateTimeArray,obj.outputs.estimatedStates(:,1))
            plot(stateTimeArray,obj.outputs.rcvrStates(:,1))

            nexttile
            hold on
            ylabel('East Velocity [m/s]')

            nexttile
            hold on
            ylabel('Down Velocity [m/s]')
            xlabel('Time [s]')


            figure
            tiledlayout(3,1)
            hold on
            ylabel('North Velocity [m/s]')

            figure
            tiledlayout(3,1)
            hold on
            ylabel('North Velocity [m/s]')
            figure
            tiledlayout(3,1)
        end

        function plotErrors(obj)
        end

        function plotCovariance(obj)
        end

        function plotResiduals(obj)
        end
    end
end

