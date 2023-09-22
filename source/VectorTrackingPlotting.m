classdef VectorTrackingPlotting < handle
    %VECTORTRACKINGPLOTTING Summary of this class goes here
    %   Detailed explanation goes here

    properties
        outputs
        timeUpdateTime
        measUpdateTime
        numMonteCarlos
        primaryColor = '#0C2340';
        secondaryColor = '#E87722';
        fs = 16;
        lw = 2;
    end

    methods
        function obj = VectorTrackingPlotting(outputFilePath)

            loadedFile = load(outputFilePath);
            obj.outputs = loadedFile.run;
            obj.timeUpdateTime = obj.outputs.estimatedStates.time;
            obj.measUpdateTime = obj.timeUpdateTime(1):1/50:obj.timeUpdateTime(end);
            obj.numMonteCarlos = 1;

            obj.plotMap();

            obj.plotStates();

            obj.plotErrors();

            obj.plotResiduals();

        end

        function plotMap(obj)
            rcvrLLA = flat2lla(obj.outputs.rcvrStates( :,7:9),[obj.outputs.trueLAT(1) obj.outputs.trueLONG(1)],0,0,'WGS84');
            estiLLA = flat2lla(obj.outputs.estimatedStates.signals.values(:,7:9),[obj.outputs.trueLAT(1) obj.outputs.trueLONG(1)],0,0,'WGS84');

            figure('Position',[200 200 900 800])
            geoplot(rcvrLLA(:,1),rcvrLLA(:,2),'LineWidth',obj.lw,'Color',obj.primaryColor)
            hold on
            geoplot(estiLLA(:,1),estiLLA(:,2),'LineWidth',obj.lw,'Color',obj.secondaryColor)
            legend('Reference','Estimated')
            ax = gca;
            ax.FontSize = obj.fs;

        end

        function plotStates(obj)
            %% NED Velocities
            figure('Position',[400 200 900 800])
            tiledlayout(3,1)
            nexttile
            hold on
            title('Velocity - NED Frame')
            ylabel('North [m/s]')
            plot(obj.measUpdateTime ,obj.outputs.rcvrStates(:,1),'LineWidth',obj.lw,'Color',obj.primaryColor)
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates.signals.values(:,1),'LineWidth',obj.lw,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('East [m/s]')
            plot(obj.measUpdateTime ,obj.outputs.rcvrStates( :,2),'LineWidth',obj.lw,'Color',obj.primaryColor)
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates.signals.values(:,2),'LineWidth',obj.lw,'Color',obj.secondaryColor)
            legend('Reference','Estimated','Location','eastoutside')
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Down [m/s]')
            xlabel('Time [s]')
            plot(obj.measUpdateTime ,obj.outputs.rcvrStates( :,3),'LineWidth',obj.lw,'Color',obj.primaryColor)
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates.signals.values(:,3),'LineWidth',obj.lw,'Color',obj.secondaryColor)
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
            plot(obj.measUpdateTime ,obj.outputs.rcvrStates( :,4)*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates.signals.values(:,4)*R2D,'LineWidth',obj.lw,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Pitch Rate [deg/s]')
            plot(obj.measUpdateTime ,obj.outputs.rcvrStates( :,5)*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates.signals.values(:,5)*R2D,'LineWidth',obj.lw,'Color',obj.secondaryColor)
            legend('Reference','Estimated','Location','eastoutside')
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Yaw Rate [deg/s]')
            xlabel('Time [s]')
            plot(obj.measUpdateTime ,obj.outputs.rcvrStates( :,6)*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates.signals.values(:,6)*R2D,'LineWidth',obj.lw,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            %% NED Positions
            figure('Position',[800 200 900 800])
            tiledlayout(3,1)
            nexttile
            hold on
            title('Position - NED Frame')
            ylabel('North [m]')
            plot(obj.measUpdateTime ,obj.outputs.rcvrStates( :,7),'LineWidth',obj.lw,'Color',obj.primaryColor)
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates.signals.values(:,7),'LineWidth',obj.lw,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('East [m]')
            plot(obj.measUpdateTime ,obj.outputs.rcvrStates( :,8),'LineWidth',obj.lw,'Color',obj.primaryColor)
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates.signals.values(:,8),'LineWidth',obj.lw,'Color',obj.secondaryColor)
            legend('Reference','Estimated','Location','eastoutside')
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Down [m]')
            xlabel('Time [s]')
            plot(obj.measUpdateTime ,obj.outputs.rcvrStates( :,9),'LineWidth',obj.lw,'Color',obj.primaryColor)
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates.signals.values(:,9),'LineWidth',obj.lw,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            %% Euler Angles
            figure('Position',[1000 200 900 800])
            tiledlayout(3,1)
            nexttile
            hold on
            title('Euler Angles')
            ylabel('Roll [deg]')
            plot(obj.measUpdateTime ,obj.outputs.rcvrStates( :,10)*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates.signals.values(:,10)*R2D,'LineWidth',obj.lw,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Pitch [deg]')
            plot(obj.measUpdateTime ,obj.outputs.rcvrStates( :,11)*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates.signals.values(:,11)*R2D,'LineWidth',obj.lw,'Color',obj.secondaryColor)
            legend('Reference','Estimated','Location','eastoutside')
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Yaw [deg]')
            xlabel('Time [s]')
            plot(obj.measUpdateTime ,obj.outputs.rcvrStates( :,12)*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates.signals.values(:,12)*R2D,'LineWidth',obj.lw,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            %% Clock Bias and Clock Drift Estimate
            figure('Position',[1200 200 900 800])
            tiledlayout(2,1)
            nexttile
            hold on
            title('Clock Bias and Drift')
            ylabel('Bias [s]')
            plot(obj.timeUpdateTime ,obj.outputs.estimatedStates.signals.values( :,13),'LineWidth',obj.lw,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Drift [s/s]')
            xlabel('Time [s]')
            plot(obj.timeUpdateTime ,obj.outputs.estimatedStates.signals.values( :,14),'LineWidth',obj.lw,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

        end

        function plotErrors(obj)
            dtMeas = diff(obj.measUpdateTime);
            dtTime = diff(obj.timeUpdateTime);

            scale = dtMeas(1)/dtTime(1);

            timeIdx = 1:scale:length(obj.timeUpdateTime);
            %% Velocity Errors
            figure('Position',[400 200 900 800])
            tiledlayout(3,1)
            nexttile
            hold on
            title('Velocity Errors')
            ylabel('North [m/s]')
            plot(obj.measUpdateTime ,obj.outputs.rcvrStates( :,1) - obj.outputs.estimatedStates.signals.values(timeIdx,1),'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('East [m/s]')
            plot(obj.measUpdateTime ,obj.outputs.rcvrStates( :,2) - obj.outputs.estimatedStates.signals.values(timeIdx,2),'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Down [m/s]')
            xlabel('Time [s]')
            plot(obj.measUpdateTime ,obj.outputs.rcvrStates( :,3) - obj.outputs.estimatedStates.signals.values(timeIdx,3),'LineWidth',obj.lw,'Color',obj.primaryColor)
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
            plot(obj.measUpdateTime ,obj.outputs.rcvrStates( :,4)*R2D - obj.outputs.estimatedStates.signals.values(timeIdx,4)*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Pitch Rate [deg/s]')
            plot(obj.measUpdateTime ,obj.outputs.rcvrStates( :,5)*R2D - obj.outputs.estimatedStates.signals.values(timeIdx,5)*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Yaw Rate [deg/s]')
            xlabel('Time [s]')
            plot(obj.measUpdateTime ,obj.outputs.rcvrStates( :,6)*R2D - obj.outputs.estimatedStates.signals.values(timeIdx,6)*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            %% Position Errors
            figure('Position',[800 200 900 800])
            tiledlayout(3,1)
            nexttile
            hold on
            title('Position Errors')
            ylabel('North [m]')
            plot(obj.measUpdateTime ,obj.outputs.rcvrStates( :,7) - obj.outputs.estimatedStates.signals.values(timeIdx,7),'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('East [m]')
            plot(obj.measUpdateTime ,obj.outputs.rcvrStates( :,8) - obj.outputs.estimatedStates.signals.values(timeIdx,8),'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Down [m]')
            xlabel('Time [s]')
            plot(obj.measUpdateTime ,obj.outputs.rcvrStates( :,9) - obj.outputs.estimatedStates.signals.values(timeIdx,9),'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            %% Euler Angle Errors
            figure('Position',[1000 200 900 800])
            tiledlayout(3,1)
            nexttile
            hold on
            title('Euler Angle Errors')
            ylabel('Roll [deg]')
            plot(obj.measUpdateTime ,obj.outputs.rcvrStates( :,10)*R2D - obj.outputs.estimatedStates.signals.values(timeIdx,10)*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Pitch [deg]')
            plot(obj.measUpdateTime ,obj.outputs.rcvrStates( :,11)*R2D - obj.outputs.estimatedStates.signals.values(timeIdx,11)*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Yaw [deg]')
            xlabel('Time [s]')
            plot(obj.measUpdateTime ,obj.outputs.rcvrStates( :,12)*R2D - obj.outputs.estimatedStates.signals.values(timeIdx,12)*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;
        end

        function plotResiduals(obj)

            avgResPsr = mean(obj.outputs.resPsr,2,"omitnan");
            avgResCarr = mean(obj.outputs.resCarr,2,"omitnan");

            figure('Position',[1400 200 900 800])
            hold on
            title('Residual Pseudoranges of In-View SVs')
            s1 = scatter(obj.measUpdateTime ,obj.outputs.resPsr,'k*');
            p1 = plot(obj.measUpdateTime ,avgResPsr ,'Color',obj.secondaryColor,LineWidth=2);
            legend([s1(1),p1(1)],'Residual Pseudoranges','Mean')
            ax = gca;
            ax.FontSize = obj.fs;

            figure('Position',[1600 200 900 800])
            hold on
            title('Residual Carrier Frequency of In-View SVs')
            s2 = scatter(obj.measUpdateTime,obj.outputs.resCarr,'k*');
            p2 = plot(obj.measUpdateTime ,avgResCarr ,'Color',obj.secondaryColor,LineWidth=2);
            legend([s2(1),p2(1)],'Residual Frequencies','Mean')
            ax = gca;
            ax.FontSize = obj.fs;

        end
    end
end

