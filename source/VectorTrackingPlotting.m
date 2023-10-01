classdef VectorTrackingPlotting < handle
    %VECTORTRACKINGPLOTTING Summary of this class goes here
    %   Detailed explanation goes here

    properties
        outputs
        timeUpdateTime
        measUpdateTime
        timeIdx
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
            obj.timeUpdateTime = obj.outputs.tout;
            obj.measUpdateTime = obj.timeUpdateTime(1):1/50:obj.timeUpdateTime(end);
            dtMeas = diff(obj.measUpdateTime);
            dtTime = diff(obj.timeUpdateTime);

            scale = dtMeas(1)/dtTime(1);

            obj.timeIdx = 1:scale:length(obj.timeUpdateTime);

            obj.numMonteCarlos = 1;

            obj.plotMap();

            obj.plotStates();

            obj.plotErrors();

            obj.plotResiduals();

            obj.plotCovariance();

            obj.plotCN0();

            obj.plotMeasurementVariance();

        end

        function plotMap(obj)
            rcvrLLA = ecef2lla(obj.outputs.rcvrStates(obj.timeIdx,1:3),'WGS84');

            figure('Position',[200 200 900 800])
            geoplot(obj.outputs.truthStates(:,7).*(180/pi),obj.outputs.truthStates(:,8).*(180/pi),'LineWidth',obj.lw,'Color','k')
            hold on
            geoplot(rcvrLLA(:,1),rcvrLLA(:,2),'LineWidth',obj.lw,'Color',obj.primaryColor)
            geoplot(obj.outputs.estimatedStates(:,7).*(180/pi),obj.outputs.estimatedStates(:,8).*(180/pi),'LineWidth',obj.lw,'Color',obj.secondaryColor)
            legend('Truth','Reference','Estimated')
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
            plot(obj.timeUpdateTime ,obj.outputs.truthStates(:,1),'LineWidth',obj.lw,'Color',obj.primaryColor)
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates(:,1),'LineWidth',obj.lw,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('East [m/s]')
            plot(obj.timeUpdateTime ,obj.outputs.truthStates(:,2),'LineWidth',obj.lw,'Color',obj.primaryColor)
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates(:,2),'LineWidth',obj.lw,'Color',obj.secondaryColor)
            legend('Truth','Estimated','Location','eastoutside')
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Down [m/s]')
            xlabel('Time [s]')
            plot(obj.timeUpdateTime ,obj.outputs.truthStates(:,3),'LineWidth',obj.lw,'Color',obj.primaryColor)
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates(:,3),'LineWidth',obj.lw,'Color',obj.secondaryColor)
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
            plot(obj.timeUpdateTime ,obj.outputs.truthStates(:,4)*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates(:,4)*R2D,'LineWidth',obj.lw,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Pitch Rate [deg/s]')
            plot(obj.timeUpdateTime ,obj.outputs.truthStates(:,5)*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates(:,5)*R2D,'LineWidth',obj.lw,'Color',obj.secondaryColor)
            legend('Truth','Estimated','Location','eastoutside')
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Yaw Rate [deg/s]')
            xlabel('Time [s]')
            plot(obj.timeUpdateTime ,obj.outputs.truthStates(:,6)*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates(:,6)*R2D,'LineWidth',obj.lw,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            %% NED Positions
            figure('Position',[800 200 900 800])
            tiledlayout(3,1)
            nexttile
            hold on
            title('Position - NED Frame')
            ylabel('North [m]')
            plot(obj.timeUpdateTime ,obj.outputs.truthStates(:,7),'LineWidth',obj.lw,'Color',obj.primaryColor)
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates(:,7),'LineWidth',obj.lw,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('East [m]')
            plot(obj.timeUpdateTime ,obj.outputs.truthStates(:,8),'LineWidth',obj.lw,'Color',obj.primaryColor)
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates(:,8),'LineWidth',obj.lw,'Color',obj.secondaryColor)
            legend('Truth','Estimated','Location','eastoutside')
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Down [m]')
            xlabel('Time [s]')
            plot(obj.timeUpdateTime ,obj.outputs.truthStates(:,9),'LineWidth',obj.lw,'Color',obj.primaryColor)
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates(:,9),'LineWidth',obj.lw,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            %% Euler Angles
            figure('Position',[1000 200 900 800])
            tiledlayout(3,1)
            nexttile
            hold on
            title('Euler Angles')
            ylabel('Roll [deg]')
            plot(obj.timeUpdateTime ,obj.outputs.truthStates(:,10)*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates(:,10)*R2D,'LineWidth',obj.lw,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Pitch [deg]')
            plot(obj.timeUpdateTime ,obj.outputs.truthStates(:,11)*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates(:,11)*R2D,'LineWidth',obj.lw,'Color',obj.secondaryColor)
            legend('Truth','Estimated','Location','eastoutside')
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Yaw [deg]')
            xlabel('Time [s]')
            plot(obj.timeUpdateTime ,obj.outputs.truthStates(:,12)*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates(:,12)*R2D,'LineWidth',obj.lw,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            %% Clock Bias and Clock Drift Estimate
            figure('Position',[1200 200 900 800])
            tiledlayout(2,1)
            nexttile
            hold on
            title('Clock Bias and Drift')
            ylabel('Bias [s]')
            plot(obj.timeUpdateTime ,obj.outputs.estimatedStates(:,13),'LineWidth',obj.lw,'Color',obj.secondaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Drift [s/s]')
            xlabel('Time [s]')
            plot(obj.timeUpdateTime ,obj.outputs.estimatedStates(:,14),'LineWidth',obj.lw,'Color',obj.secondaryColor)
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
            plot(obj.timeUpdateTime ,obj.outputs.truthStates(:,1) - obj.outputs.estimatedStates(:,1),'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('East [m/s]')
            plot(obj.timeUpdateTime ,obj.outputs.truthStates(:,2) - obj.outputs.estimatedStates(:,2),'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Down [m/s]')
            xlabel('Time [s]')
            plot(obj.timeUpdateTime ,obj.outputs.truthStates(:,3) - obj.outputs.estimatedStates(:,3),'LineWidth',obj.lw,'Color',obj.primaryColor)
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
            plot(obj.timeUpdateTime ,(obj.outputs.truthStates(:,4) - obj.outputs.estimatedStates(:,4))*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Pitch Rate [deg/s]')
            plot(obj.timeUpdateTime ,(obj.outputs.truthStates(:,5) - obj.outputs.estimatedStates(:,5))*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Yaw Rate [deg/s]')
            xlabel('Time [s]')
            plot(obj.timeUpdateTime ,(obj.outputs.truthStates(:,6) - obj.outputs.estimatedStates(:,6))*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            %% Position Errors
            figure('Position',[800 200 900 800])
            tiledlayout(3,1)
            nexttile
            hold on
            title('Position Errors')
            ylabel('North [m]')
            plot(obj.timeUpdateTime ,(obj.outputs.truthStates(:,7) - obj.outputs.estimatedStates(:,7))*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('East [m]')
            plot(obj.timeUpdateTime ,(obj.outputs.truthStates(:,8) - obj.outputs.estimatedStates(:,8))*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Down [m]')
            xlabel('Time [s]')
            plot(obj.timeUpdateTime ,obj.outputs.truthStates(:,9) - obj.outputs.estimatedStates(:,9),'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            %% Euler Angle Errors
            figure('Position',[1000 200 900 800])
            tiledlayout(3,1)
            nexttile
            hold on
            title('Euler Angle Errors')
            ylabel('Roll [deg]')
            plot(obj.timeUpdateTime ,(obj.outputs.truthStates(:,10) - obj.outputs.estimatedStates(:,10))*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Pitch [deg]')
            plot(obj.timeUpdateTime ,(obj.outputs.truthStates(:,11) - obj.outputs.estimatedStates(:,11))*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Yaw [deg]')
            xlabel('Time [s]')
            plot(obj.timeUpdateTime ,(obj.outputs.truthStates(:,12) - obj.outputs.estimatedStates(:,12))*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;
        end

        function plotResiduals(obj)

            avgResPsr = mean(obj.outputs.resPsr,2,"omitnan");
            avgResCarr = mean(obj.outputs.resCarr,2,"omitnan");

            figure('Position',[1400 200 900 800])
            hold on
            title('Residual Pseudoranges of In-View SVs')
            s1 = scatter(obj.measUpdateTime ,obj.outputs.resPsr(obj.timeIdx,:),'k*');
            p1 = plot(obj.measUpdateTime ,avgResPsr(obj.timeIdx) ,'Color',obj.secondaryColor,LineWidth=2);
            legend([s1(1),p1(1)],'Residual Pseudoranges','Mean')
            ax = gca;
            ax.FontSize = obj.fs;

            figure('Position',[1600 200 900 800])
            hold on
            title('Residual Carrier Frequency of In-View SVs')
            s2 = scatter(obj.measUpdateTime,obj.outputs.resCarr(obj.timeIdx,:),'k*');
            p2 = plot(obj.measUpdateTime ,avgResCarr(obj.timeIdx) ,'Color',obj.secondaryColor,LineWidth=2);
            legend([s2(1),p2(1)],'Residual Frequencies','Mean')
            ax = gca;
            ax.FontSize = obj.fs;

        end

        function plotCovariance(obj)
            for i = 1:size(obj.outputs.covariance,3)
                velocity(i,:) = trace(obj.outputs.covariance(1:3,1:3,i));
                omega(i,:) = trace(obj.outputs.covariance(4:6,4:6,i));
                position(i,:) = trace(obj.outputs.covariance(7:9,7:9,i));
                angles(i,:) = trace(obj.outputs.covariance(10:12,10:12,i));
                clkBias(i,:) = (obj.outputs.covariance(13,13,i));
                clkDrift(i,:) = (obj.outputs.covariance(14,14,i));

            end
            figure('Position',[400 0 900 1600])
            tiledlayout(5,1)
            nexttile
            hold on
            title('Covariance Magnitudes')
            ylabel('Velocity [m/s]')
            plot(obj.timeUpdateTime ,velocity,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Angular Rate [rad/s]')
            plot(obj.timeUpdateTime ,omega,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Position [m]')
            plot(obj.timeUpdateTime ,position,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Euler Angles [rad]')
            plot(obj.timeUpdateTime ,angles,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Clock Terms')
            xlabel('Time [s]')
            plot(obj.timeUpdateTime ,clkBias,'LineWidth',obj.lw,'Color',obj.primaryColor)
            plot(obj.timeUpdateTime ,clkDrift,'LineWidth',obj.lw,'Color',obj.secondaryColor)
            legend('Bias [s]','Drift [s/s]')
            ax = gca;
            ax.FontSize = obj.fs;

        end

        function plotCN0(obj)

            CN0idx = diff(obj.outputs.CN0,1,1) ~= 0;
            CN0idx = find(CN0idx(obj.timeIdx(2)-1,:));
            for i = 1:length(CN0idx)
            CN0plot(:,i) = 10*log(obj.outputs.CN0(:,CN0idx(i)))/log(10);

            end

            if ~isempty(CN0idx)
                figure('Position',[200 0 900 800])
                hold on
                title('Carrier to Noise Ratio')
                ylabel('CN0 [db-Hz]')

                plot(obj.timeUpdateTime(obj.timeIdx),CN0plot(obj.timeIdx,:),'LineWidth',obj.lw,'Color',obj.primaryColor)
                ax = gca;
                ax.FontSize = obj.fs;
            end

        end

        function plotMeasurementVariance(obj)

            avgVarPsr = mean(obj.outputs.varPsr,2,"omitnan");
            avgVarCarr = mean(obj.outputs.varCarr,2,"omitnan");

            figure('Position',[1400 200 900 800])
            hold on
            title('Pseudorange Variance of In-View SVs')
            s1 = scatter(obj.measUpdateTime ,obj.outputs.varPsr(obj.timeIdx,:),'k*');
            p1 = plot(obj.measUpdateTime ,avgVarPsr(obj.timeIdx) ,'Color',obj.secondaryColor,LineWidth=2);
            legend([s1(1),p1(1)],'Pseudorange Variance','Mean')
            ylim([0 500])
            ax = gca;
            ax.FontSize = obj.fs;

            figure('Position',[1600 200 900 800])
            hold on
            title('Pseudorange-Rate Variance of In-View SVs')
            s2 = scatter(obj.measUpdateTime,obj.outputs.varCarr(obj.timeIdx,:),'k*');
            p2 = plot(obj.measUpdateTime ,avgVarCarr(obj.timeIdx) ,'Color',obj.secondaryColor,LineWidth=2);
            legend([s2(1),p2(1)],'Pseudorange-Rate Variance','Mean')
            ylim([0 500])
            ax = gca;
            ax.FontSize = obj.fs;
        end
    end
end

