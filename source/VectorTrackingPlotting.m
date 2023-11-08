classdef VectorTrackingPlotting < handle

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
            obj.outputs = loadedFile.run(1,1);
            obj.timeUpdateTime = obj.outputs.tout;
            obj.measUpdateTime = obj.timeUpdateTime(1):1/50:obj.timeUpdateTime(end);
            dtMeas = diff(obj.measUpdateTime);
            dtTime = diff(obj.timeUpdateTime);

            scale = dtMeas(1)/dtTime(1);

            obj.timeIdx = 1:scale:length(obj.timeUpdateTime);

            obj.plotMap();

            obj.plotStates();

            % obj.plotErrors();

            obj.plotResiduals();

            obj.plotCN0();

            obj.plotMeasurementVariance();

            obj.satellitePlot();

            % for i = [1 7 13 14 15 17 19 21 29]
                figure('Position',[200 200 900 400])
                hold on
                s1 = plot(movmean(sqrt(obj.outputs.IE(1:4:end,[1 7 13 14 15 17 19 21 29]).^2 + obj.outputs.QE(1:4:end,[1 7 13 14 15 17 19 21 29]).^2),100),'Color',obj.primaryColor,'LineWidth',2);
                s2 = plot((movmean(sqrt(obj.outputs.IP(1:4:end,[1 7 13 14 15 17 19 21 29]).^2 + obj.outputs.QP(1:4:end,[1 7 13 14 15 17 19 21 29]).^2),100)),'Color',obj.secondaryColor,'LineWidth',2);
                s3 = plot(movmean(sqrt(obj.outputs.IL(1:4:end,[1 7 13 14 15 17 19 21 29]).^2 + obj.outputs.QL(1:4:end,[1 7 13 14 15 17 19 21 29]).^2),100),'Color','k','LineWidth',2);
                legend([s1(1) s2(1) s3(1)],'Early Power','Prompt Power','Late Power')
                xlabel('Samples')
                xlim([0 3000])
                ax = gca;
                ax.FontSize = obj.fs;
            % end

            for i = [1 7 13 14 15 17 19 21 29]
            figure('Position',[200 200 900 800])
            hold on

            s2 = plot(obj.outputs.IP(1:4:end,i),obj.outputs.QP(1:4:end,i),'*');
            end
            residualCodePhase = obj.outputs.resPsr(1:4:end,[1 7 13 14 15 17 19 21 29])/(299792458/1.023e6);
            residualCarrierFreq = obj.outputs.resCarr(1:4:end,[1 7 13 14 15 17 19 21 29])/(-299792458/1575.42e6);

            disp(mean(rms(residualCodePhase)))
            disp(mean(rms(residualCarrierFreq)))
        end

        function plotMap(obj)

            figure('Position',[200 200 900 800])
            geoplot(obj.outputs.estimatedStates(:,7).*(180/pi),obj.outputs.estimatedStates(:,8).*(180/pi),'*','Color',obj.secondaryColor)
            hold on
            geoplot(obj.outputs.truthStates(:,7).*(180/pi),obj.outputs.truthStates(:,8).*(180/pi),'LineWidth',obj.lw,'Color',obj.primaryColor)
            legend('Estimated','Truth')
            ax = gca;
            ax.FontSize = 20;

        end

        function plotStates(obj)
            %% NED Velocities
            figure('Position',[400 200 900 800])
            tiledlayout(3,1)
            nexttile
            hold on
            title('Velocity - NED Frame')
            ylabel('North [m/s]')
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates(:,1),'LineWidth',obj.lw,'Color',obj.secondaryColor)
            plot(obj.measUpdateTime ,obj.outputs.truthStates(:,1),'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('East [m/s]')
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates(:,2),'LineWidth',obj.lw,'Color',obj.secondaryColor)
            plot(obj.measUpdateTime ,obj.outputs.truthStates(:,2),'LineWidth',obj.lw,'Color',obj.primaryColor)
            legend('Estimated','Truth','Location','eastoutside')
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Down [m/s]')
            xlabel('Time [s]')
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates(:,3),'LineWidth',obj.lw,'Color',obj.secondaryColor)
            plot(obj.measUpdateTime ,obj.outputs.truthStates(:,3),'LineWidth',obj.lw,'Color',obj.primaryColor)
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
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates(:,4)*R2D,'LineWidth',obj.lw,'Color',obj.secondaryColor)
            plot(obj.measUpdateTime ,obj.outputs.truthStates(:,4)*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ylim([-180 180])
            yticks([-180 -90 0 90 180])
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Pitch Rate [deg/s]')
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates(:,5)*R2D,'LineWidth',obj.lw,'Color',obj.secondaryColor)
            plot(obj.measUpdateTime ,obj.outputs.truthStates(:,5)*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ylim([-180 180])
            yticks([-180 -90 0 90 180])
            legend('Estimated','Truth','Location','eastoutside')
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Yaw Rate [deg/s]')
            xlabel('Time [s]')
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates(:,6)*R2D,'LineWidth',obj.lw,'Color',obj.secondaryColor)
            plot(obj.measUpdateTime ,obj.outputs.truthStates(:,6)*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ylim([-180 180])
            yticks([-180 -90 0 90 180])
            ax = gca;
            ax.FontSize = obj.fs;

            %% NED Positions
            truthNED = lla2ned(obj.outputs.truthStates(:,7:9),[obj.outputs.truthStates(1,7:8), 0],"ellipsoid");
            noiseNED = lla2ned(obj.outputs.estimatedStates(:,7:9),[obj.outputs.estimatedStates(1,7:8), 0],"ellipsoid");
            figure('Position',[800 200 900 800])
            tiledlayout(3,1)
            nexttile
            hold on
            title('Position - NED Frame')
            ylabel('North [m]')
            plot(obj.timeUpdateTime,noiseNED(:,1),'LineWidth',obj.lw,'Color',obj.secondaryColor)
            plot(obj.measUpdateTime ,truthNED(:,1),'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('East [m]')
            plot(obj.timeUpdateTime,noiseNED(:,2),'LineWidth',obj.lw,'Color',obj.secondaryColor)
            plot(obj.measUpdateTime ,truthNED(:,2),'LineWidth',obj.lw,'Color',obj.primaryColor)
            legend('Estimated','Truth','Location','eastoutside')
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Down [m]')
            xlabel('Time [s]')
            plot(obj.timeUpdateTime,noiseNED(:,3),'LineWidth',obj.lw,'Color',obj.secondaryColor)
            plot(obj.measUpdateTime ,truthNED(:,3),'LineWidth',obj.lw,'Color',obj.primaryColor)
            xlim([0 180])
            ax = gca;
            ax.FontSize = obj.fs;

            %% Euler Angles
            figure('Position',[1000 200 900 800])
            tiledlayout(3,1)
            nexttile
            hold on
            title('Euler Angles')
            ylabel('Roll [deg]')
            plot(obj.timeUpdateTime,wrapTo180(obj.outputs.estimatedStates(:,10)*R2D),'LineWidth',obj.lw,'Color',obj.secondaryColor)
            plot(obj.measUpdateTime ,obj.outputs.truthStates(:,10)*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ylim([-180 180])
            yticks([-180 -90 0 90 180])
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Pitch [deg]')
            plot(obj.timeUpdateTime,wrapTo180(obj.outputs.estimatedStates(:,11)*R2D),'LineWidth',obj.lw,'Color',obj.secondaryColor)
            plot(obj.measUpdateTime ,obj.outputs.truthStates(:,11)*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            legend('Estimated','Truth','Location','eastoutside')
            ylim([-180 180])
            yticks([-180 -90 0 90 180])
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Yaw [deg]')
            xlabel('Time [s]')
            plot(obj.timeUpdateTime,wrapTo180(obj.outputs.estimatedStates(:,12)*R2D),'LineWidth',obj.lw,'Color',obj.secondaryColor)
            plot(obj.measUpdateTime ,obj.outputs.truthStates(:,12)*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ylim([-180 180])
            yticks([-180 -90 0 90 180])
            ax = gca;
            ax.FontSize = obj.fs;

            %% Clock Bias and Clock Drift Estimate
            figure('Position',[1200 200 900 800])
            tiledlayout(2,1)
            nexttile
            hold on
            title('Clock Bias and Drift')
            ylabel('Bias [m]')
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates(:,13),'LineWidth',obj.lw,'Color',obj.secondaryColor)
            plot(obj.measUpdateTime ,obj.outputs.truthStates(:,13),'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Drift [m/s]')
            xlabel('Time [s]')
            plot(obj.timeUpdateTime,obj.outputs.estimatedStates(:,14),'LineWidth',obj.lw,'Color',obj.secondaryColor)
            plot(obj.measUpdateTime ,obj.outputs.truthStates(:,14),'LineWidth',obj.lw,'Color',obj.primaryColor)

            legend('Estimated','Truth','Location','northoutside')
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
            ylim([-180 180])
            yticks([-180 -90 0 90 180])
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Pitch Rate [deg/s]')
            plot(obj.timeUpdateTime ,(obj.outputs.truthStates(:,5) - obj.outputs.estimatedStates(:,5))*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ylim([-180 180])
            yticks([-180 -90 0 90 180])
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Yaw Rate [deg/s]')
            xlabel('Time [s]')
            plot(obj.timeUpdateTime ,(obj.outputs.truthStates(:,6) - obj.outputs.estimatedStates(:,6))*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ylim([-180 180])
            yticks([-180 -90 0 90 180])
            ax = gca;
            ax.FontSize = obj.fs;

            %% Position Errors
            truthNED = lla2ned(obj.outputs.truthStates(:,7:9),obj.outputs.truthStates(1,7:9),"flat");
            noiseNED = lla2ned(obj.outputs.estimatedStates(:,7:9),obj.outputs.estimatedStates(1,7:9),"flat");
            figure('Position',[800 200 900 800])
            tiledlayout(3,1)
            nexttile
            hold on
            title('Position Errors')
            ylabel('North [m]')
            plot(obj.timeUpdateTime ,(truthNED(:,1) - noiseNED(:,1)),'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('East [m]')
            plot(obj.timeUpdateTime ,(truthNED(:,2) - noiseNED(:,2)),'LineWidth',obj.lw,'Color',obj.primaryColor)
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Down [m]')
            xlabel('Time [s]')
            plot(obj.timeUpdateTime ,(truthNED(:,3) - noiseNED(:,3)),'LineWidth',obj.lw,'Color',obj.primaryColor)
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
            ylim([-180 180])
            yticks([-180 -90 0 90 180])
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Pitch [deg]')
            plot(obj.timeUpdateTime ,(obj.outputs.truthStates(:,11) - obj.outputs.estimatedStates(:,11))*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ylim([-180 180])
            yticks([-180 -90 0 90 180])
            ax = gca;
            ax.FontSize = obj.fs;

            nexttile
            hold on
            ylabel('Yaw [deg]')
            xlabel('Time [s]')
            plot(obj.timeUpdateTime ,(obj.outputs.truthStates(:,12) - obj.outputs.estimatedStates(:,12))*R2D,'LineWidth',obj.lw,'Color',obj.primaryColor)
            ylim([-180 180])
            yticks([-180 -90 0 90 180])
            ax = gca;
            ax.FontSize = obj.fs;
        end

        function plotResiduals(obj)

            avgResPsr = mean(obj.outputs.resPsr,2,"omitnan");
            avgResCarr = mean(obj.outputs.resCarr,2,"omitnan");

            if ~isnan(avgResCarr(1))
                figure('Position',[1400 200 900 800])
                hold on
                % title('Residual Pseudoranges of In-View SVs')
                s1 = scatter(obj.measUpdateTime ,obj.outputs.resPsr(obj.timeIdx,:),'k*');
                p1 = plot(obj.measUpdateTime ,avgResPsr(obj.timeIdx) ,'Color',obj.secondaryColor,LineWidth=2);
                xlabel('Time [s]')
                ylabel('Residuals [m]')
                legend([s1(1),p1(1)],'Residual Pseudoranges','Mean')
                % ylim([0 100])
                xlim([60 180])
                ax = gca;
                ax.FontSize = obj.fs;

                figure('Position',[1600 200 900 800])
                hold on
                % title('Residual Pseudorange-Rate of In-View SVs')
                s2 = scatter(obj.measUpdateTime,obj.outputs.resCarr(obj.timeIdx,:),'k*');
                p2 = plot(obj.measUpdateTime ,avgResCarr(obj.timeIdx) ,'Color',obj.secondaryColor,LineWidth=2);
                xlabel('Time [s]')
                ylabel('Residuals [m/s]')
                legend([s2(1),p2(1)],'Residual Pseudorange-Rates','Mean')
                % ylim([0 100])
                xlim([60 180])
                ax = gca;
                ax.FontSize = obj.fs;
            else
                printText(12);
            end

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

            % CN0idx = diff(obj.outputs.CN0,1,1) ~= 0;
            CN0idx = [1 7 13 14 15 17 19 21 29];
            for i = 1:length(CN0idx)
                CN0plot(:,i) = 10*log(obj.outputs.CN0(:,CN0idx(i)))/log(10);

            end

            if ~isempty(CN0idx)
                figure('Position',[200 0 900 800])
                hold on
                title('Carrier to Noise Ratio')
                ylabel('CN0 [db-Hz]')
                xlabel('Time [s]')
                first = 0:1/200:60;
                second = 60+1/200:1/200:150;
                third = 150+1/200:1/200:180;

                s1 = plot(obj.timeUpdateTime(obj.timeIdx),CN0plot(obj.timeIdx,:),'LineWidth',obj.lw,'Color',obj.primaryColor);
                % s2 = plot(obj.timeUpdateTime,[ones(1,length(first)).*45 ones(1,length(second)).*25 ones(1,length(third)).*45],'LineWidth',3,'Color',obj.secondaryColor,'LineStyle','--');
                % legend([s1(1),s2(1)],'Estimated CN0','Design CN0')
                ylim([0 50])
                ax = gca;
                ax.FontSize = obj.fs;
            else
                printText(12);
            end

        end

        function plotMeasurementVariance(obj)


            avgVarPsr = mean(obj.outputs.varPsr,2,"omitnan");
            avgVarCarr = mean(obj.outputs.varCarr,2,"omitnan");
            minVarPsr = min(min(obj.outputs.varPsr));
            minVarCarr = min(min(obj.outputs.varCarr));
            maxVarPsr = max(max(obj.outputs.varPsr(1000:end,:)));
            maxVarCarr = max(max(obj.outputs.varCarr(1000:end,:)));
            if ~isnan(avgVarPsr(1))

                figure('Position',[1400 200 900 800])
                hold on
                % title('Pseudorange Variance of In-View SVs')
                s1 = scatter(obj.measUpdateTime ,obj.outputs.varPsr(obj.timeIdx,:),'k*');
                p1 = plot(obj.measUpdateTime ,avgVarPsr(obj.timeIdx) ,'Color',obj.secondaryColor,LineWidth=2);
                legend([s1(1),p1(1)],'Pseudorange Variance','Mean')
                xlabel('Time [s]')
                ylabel('Variance [m]')
                % ylim([0 10000])
                xlim([60 180])
                ax = gca;
                ax.FontSize = obj.fs;

                figure('Position',[1600 200 900 800])
                hold on
                % title('Pseudorange-Rate Variance of In-View SVs')
                s2 = scatter(obj.measUpdateTime,obj.outputs.varCarr(obj.timeIdx,:),'k*');
                p2 = plot(obj.measUpdateTime ,avgVarCarr(obj.timeIdx) ,'Color',obj.secondaryColor,LineWidth=2);
                legend([s2(1),p2(1)],'Pseudorange-Rate Variance','Mean')
                xlabel('Time [s]')
                ylabel('Variance [m/s]')
                % ylim([0 100])
                xlim([60 180])
                ax = gca;
                ax.FontSize = obj.fs;
            else
                printText(12);
            end
        end

        function satellitePlot(obj)
            figure('Position',[1400 200 900 800])

            truthLat = obj.outputs.truthStates(:,7).*(180/pi);
            truthLong = obj.outputs.truthStates(:,8).*(180/pi);
            truthalt = obj.outputs.truthStates(:,9);

            noiseLat = obj.outputs.estimatedStates(:,7).*(180/pi);
            noiseLong = obj.outputs.estimatedStates(:,8).*(180/pi);
            noisealt = obj.outputs.estimatedStates(:,9);
            count = 1;
            for i = obj.timeIdx
                [svaz,svel,~] = ecef2aer(obj.outputs.svStates(1,:,i),obj.outputs.svStates(3,:,i),obj.outputs.svStates(5,:,i),truthLat(i),truthLong(i),truthalt(i),wgs84Ellipsoid("meter"),"degrees");
                count = count + 1;
                idx = svel > 10;
                svel_truth(:,count) = svel(idx);
                svaz_truth(:,count) = svaz(idx);

            end

            sp = skyplot(svaz_truth(:,2),svel_truth(:,2));

        end
    end
end

