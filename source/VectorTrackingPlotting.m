classdef VectorTrackingPlotting < handle
    %VECTORTRACKINGPLOTTING Summary of this class goes here
    %   Detailed explanation goes here

    properties
        outputs
        timeVector
        numMonteCarlos
        primaryColor = '#0C2340';
        secondaryColor = '#E87722';
        fs = 16;
    end

    methods
        function obj = VectorTrackingPlotting(outputFilePath)

            loadedFile = load(outputFilePath);
            obj.outputs = loadedFile.run;
            obj.timeVector = obj.outputs.estimatedStates.time;
            obj.numMonteCarlos = 1;

            obj.plotStates();

            obj.plotErrors();

            obj.plotResiduals();




        end

        function plotStates(obj)
            %% NED Velocities
            figure('Position',[400 200 900 800])
            tiledlayout(3,1)
            nexttile
            hold on
            title('Velocity - NED Frame')
            ylabel('North [m/s]')
            plot(obj.timeVector(1:4:end),obj.outputs.rcvrStates(1:4:end,1),'LineWidth',1,'Color',obj.primaryColor)
            plot(obj.timeVector,obj.outputs.estimatedStates.signals.values(:,1),'LineWidth',1,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('East [m/s]')
            plot(obj.timeVector(1:4:end),obj.outputs.rcvrStates(1:4:end,2),'LineWidth',1,'Color',obj.primaryColor)
            plot(obj.timeVector,obj.outputs.estimatedStates.signals.values(:,2),'LineWidth',1,'Color',obj.secondaryColor)
            legend('Reference','Estimated','Location','eastoutside')
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Down [m/s]')
            xlabel('Time [s]')
            plot(obj.timeVector(1:4:end),obj.outputs.rcvrStates(1:4:end,3),'LineWidth',1,'Color',obj.primaryColor)
            plot(obj.timeVector,obj.outputs.estimatedStates.signals.values(:,3),'LineWidth',1,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            %% Angular Rates
            R2D = (180/pi);
            figure('Position',[600 200 900 800])
            tiledlayout(3,1)
            nexttile
            hold on
            title('Angular Rates - Body Frame')
            ylabel('Roll Rate [deg/s]')
            plot(obj.timeVector(1:4:end),obj.outputs.rcvrStates(1:4:end,4)*R2D,'LineWidth',1,'Color',obj.primaryColor)
            plot(obj.timeVector,obj.outputs.estimatedStates.signals.values(:,4)*R2D,'LineWidth',1,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Pitch Rate [deg/s]')
            plot(obj.timeVector(1:4:end),obj.outputs.rcvrStates(1:4:end,5)*R2D,'LineWidth',1,'Color',obj.primaryColor)
            plot(obj.timeVector,obj.outputs.estimatedStates.signals.values(:,5)*R2D,'LineWidth',1,'Color',obj.secondaryColor)
            legend('Reference','Estimated','Location','eastoutside')
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Yaw Rate [deg/s]')
            xlabel('Time [s]')
            plot(obj.timeVector(1:4:end),obj.outputs.rcvrStates(1:4:end,6)*R2D,'LineWidth',1,'Color',obj.primaryColor)
            plot(obj.timeVector,obj.outputs.estimatedStates.signals.values(:,6)*R2D,'LineWidth',1,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            %% NED Positions
            figure('Position',[800 200 900 800])
            tiledlayout(3,1)
            nexttile
            hold on
            title('Position - NED Frame')
            ylabel('North [m]')
            plot(obj.timeVector(1:4:end),obj.outputs.rcvrStates(1:4:end,7),'LineWidth',1,'Color',obj.primaryColor)
            plot(obj.timeVector,obj.outputs.estimatedStates.signals.values(:,7),'LineWidth',1,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('East [m]')
            plot(obj.timeVector(1:4:end),obj.outputs.rcvrStates(1:4:end,8),'LineWidth',1,'Color',obj.primaryColor)
            plot(obj.timeVector,obj.outputs.estimatedStates.signals.values(:,8),'LineWidth',1,'Color',obj.secondaryColor)
            legend('Reference','Estimated','Location','eastoutside')
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Down [m]')
            xlabel('Time [s]')
            plot(obj.timeVector(1:4:end),obj.outputs.rcvrStates(1:4:end,9),'LineWidth',1,'Color',obj.primaryColor)
            plot(obj.timeVector,obj.outputs.estimatedStates.signals.values(:,9),'LineWidth',1,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            %% Euler Angles
            figure('Position',[1000 200 900 800])
            tiledlayout(3,1)
            nexttile
            hold on
            title('Euler Angles')
            ylabel('Roll [deg]')
            plot(obj.timeVector(1:4:end),obj.outputs.rcvrStates(1:4:end,10)*R2D,'LineWidth',1,'Color',obj.primaryColor)
            plot(obj.timeVector,obj.outputs.estimatedStates.signals.values(:,10)*R2D,'LineWidth',1,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Pitch [deg]')
            plot(obj.timeVector(1:4:end),obj.outputs.rcvrStates(1:4:end,11)*R2D,'LineWidth',1,'Color',obj.primaryColor)
            plot(obj.timeVector,obj.outputs.estimatedStates.signals.values(:,11)*R2D,'LineWidth',1,'Color',obj.secondaryColor)
            legend('Reference','Estimated','Location','eastoutside')
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Yaw [deg]')
            xlabel('Time [s]')
            plot(obj.timeVector(1:4:end),obj.outputs.rcvrStates(1:4:end,12)*R2D,'LineWidth',1,'Color',obj.primaryColor)
            plot(obj.timeVector,obj.outputs.estimatedStates.signals.values(:,12)*R2D,'LineWidth',1,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            %% Clock Bias and Clock Drift Estimate
            figure('Position',[1200 200 900 800])
            tiledlayout(2,1)
            nexttile
            hold on
            title('Clock Bias and Drift')
            ylabel('Bias [s]')
            plot(obj.timeVector,obj.outputs.estimatedStates.signals.values(:,13),'LineWidth',1,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Drift [s/s]')
            xlabel('Time [s]')
            plot(obj.timeVector,obj.outputs.estimatedStates.signals.values(:,14),'LineWidth',1,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

        end

        function plotErrors(obj)

            %% Velocity Errors
            figure('Position',[400 200 900 800])
            tiledlayout(3,1)
            nexttile
            hold on
            title('Velocity Errors')
            ylabel('North [m/s]')
            plot(obj.timeVector(1:4:end),obj.outputs.rcvrStates(1:4:end,1) - obj.outputs.estimatedStates.signals.values(1:4:end,1),'LineWidth',1,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('East [m/s]')
            plot(obj.timeVector(1:4:end),obj.outputs.rcvrStates(1:4:end,2) - obj.outputs.estimatedStates.signals.values(1:4:end,2),'LineWidth',1,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Down [m/s]')
            xlabel('Time [s]')
            plot(obj.timeVector(1:4:end),obj.outputs.rcvrStates(1:4:end,3) - obj.outputs.estimatedStates.signals.values(1:4:end,3),'LineWidth',1,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            %% Angular Rate Errors
            R2D = (180/pi);
            figure('Position',[600 200 900 800])
            tiledlayout(3,1)
            nexttile
            hold on
            title('Angular Rate Errors')
            ylabel('Roll Rate [deg/s]')
            plot(obj.timeVector(1:4:end),obj.outputs.rcvrStates(1:4:end,4)*R2D - obj.outputs.estimatedStates.signals.values(1:4:end,4)*R2D,'LineWidth',1,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Pitch Rate [deg/s]')
            plot(obj.timeVector(1:4:end),obj.outputs.rcvrStates(1:4:end,5)*R2D - obj.outputs.estimatedStates.signals.values(1:4:end,5)*R2D,'LineWidth',1,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Yaw Rate [deg/s]')
            xlabel('Time [s]')
            plot(obj.timeVector(1:4:end),obj.outputs.rcvrStates(1:4:end,6)*R2D - obj.outputs.estimatedStates.signals.values(1:4:end,6)*R2D,'LineWidth',1,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            %% Position Errors
            figure('Position',[800 200 900 800])
            tiledlayout(3,1)
            nexttile
            hold on
            title('Position Errors')
            ylabel('North [m]')
            plot(obj.timeVector(1:4:end),obj.outputs.rcvrStates(1:4:end,7) - obj.outputs.estimatedStates.signals.values(1:4:end,7),'LineWidth',1,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('East [m]')
            plot(obj.timeVector(1:4:end),obj.outputs.rcvrStates(1:4:end,8) - obj.outputs.estimatedStates.signals.values(1:4:end,8),'LineWidth',1,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Down [m]')
            xlabel('Time [s]')
            plot(obj.timeVector(1:4:end),obj.outputs.rcvrStates(1:4:end,9) - obj.outputs.estimatedStates.signals.values(1:4:end,9),'LineWidth',1,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            %% Euler Angle Errors
            figure('Position',[1000 200 900 800])
            tiledlayout(3,1)
            nexttile
            hold on
            title('Euler Angle Errors')
            ylabel('Roll [deg]')
            plot(obj.timeVector(1:4:end),obj.outputs.rcvrStates(1:4:end,10)*R2D - obj.outputs.estimatedStates.signals.values(1:4:end,10)*R2D,'LineWidth',1,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Pitch [deg]')
            plot(obj.timeVector(1:4:end),obj.outputs.rcvrStates(1:4:end,11)*R2D - obj.outputs.estimatedStates.signals.values(1:4:end,11)*R2D,'LineWidth',1,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Yaw [deg]')
            xlabel('Time [s]')
            plot(obj.timeVector(1:4:end),obj.outputs.rcvrStates(1:4:end,12)*R2D - obj.outputs.estimatedStates.signals.values(1:4:end,12)*R2D,'LineWidth',1,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;
        end

        function plotResiduals(obj)

            avgResPsr = mean(obj.outputs.resPsr,2,"omitmissing");
            avgResCarr = mean(obj.outputs.resCarr,2,"omitmissing");

            figure('Position',[1400 200 900 800])
            hold on
            title('Residual Pseudoranges of In-View SVs')
            s1 = scatter(obj.timeVector,obj.outputs.resPsr,'k*');
            p1 = plot(obj.timeVector(1:4:end),avgResPsr(1:4:end),'Color',obj.secondaryColor,LineWidth=2);
            legend([s1(1),p1(1)],'Residual Pseudoranges','Mean')
            ax = gca;
            ax.FontSize = obj.fs;

            figure('Position',[1600 200 900 800])
            hold on
            title('Residual Carrier Frequency of In-View SVs')
            s2 = scatter(obj.timeVector,obj.outputs.resCarr,'k*');
            p2 = plot(obj.timeVector(1:4:end),avgResCarr(1:4:end),'Color',obj.secondaryColor,LineWidth=2);
            legend([s2(1),p2(1)],'Residual Frequencies','Mean')
            ax = gca;
            ax.FontSize = obj.fs;

        end
    end
end

