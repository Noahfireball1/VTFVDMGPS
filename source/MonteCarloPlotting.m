classdef MonteCarloPlotting < handle

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
        function obj = MonteCarloPlotting(outputFilePath)
            mcGoodIdx = [];
            loadedFile = load(outputFilePath);
            numMCRuns = size(loadedFile.run,2);
            % Determining RMS value for each state
            for i = 1:100
                successfulRuns = loadedFile.run(1,i).estimatedStates(end,9);

                if ~isnan(successfulRuns) && successfulRuns > 0
                    mcGoodIdx = [mcGoodIdx i];
                end
            end


            count = 1;
            for timeIdx = 1:20:length(loadedFile.run(1,1).tout)
                for mcIdx = [1:7 9:100]
                    u(count,mcIdx) = loadedFile.run(1,mcIdx).estimatedStates(timeIdx,1);
                    v(count,mcIdx)= loadedFile.run(1,mcIdx).estimatedStates(timeIdx,2);
                    w(count,mcIdx) = loadedFile.run(1,mcIdx).estimatedStates(timeIdx,3);
                    speed(count,mcIdx)= norm([u(count) v(count) w(count)]);
                    p(count,mcIdx) = loadedFile.run(1,mcIdx).estimatedStates(timeIdx,4);
                    q(count,mcIdx) = loadedFile.run(1,mcIdx).estimatedStates(timeIdx,5);
                    r(count,mcIdx)= loadedFile.run(1,mcIdx).estimatedStates(timeIdx,6);
                    lat(count,mcIdx) = loadedFile.run(1,mcIdx).estimatedStates(timeIdx,7);
                    long(count,mcIdx) = loadedFile.run(1,mcIdx).estimatedStates(timeIdx,8);
                    alt(count,mcIdx) = loadedFile.run(1,mcIdx).estimatedStates(timeIdx,9);
                    phi(count,mcIdx) = loadedFile.run(1,mcIdx).estimatedStates(timeIdx,10);
                    theta(count,mcIdx) = loadedFile.run(1,mcIdx).estimatedStates(timeIdx,11);
                    psi(count,mcIdx) = loadedFile.run(1,mcIdx).estimatedStates(timeIdx,12);
                    clkBias(count,mcIdx) = loadedFile.run(1,mcIdx).estimatedStates(timeIdx,13);
                    clkDrift(count,mcIdx) = loadedFile.run(1,mcIdx).estimatedStates(timeIdx,14);
                    codephaseError(mcIdx,:,count) = loadedFile.run(1,mcIdx).resPsr(timeIdx,:);
                    carrFreqError(mcIdx,:,count) = loadedFile.run(1,mcIdx).resCarr(timeIdx,:);

                    true_u(count,mcIdx) = loadedFile.run(1,mcIdx).truthStates(timeIdx,1);
                    true_v(count,mcIdx) = loadedFile.run(1,mcIdx).truthStates(timeIdx,2);
                    true_w(count,mcIdx) = loadedFile.run(1,mcIdx).truthStates(timeIdx,3);
                    true_speed(count,mcIdx) = norm([true_u(count) true_v(count) true_w(count)]);
                    true_p(count,mcIdx) = loadedFile.run(1,mcIdx).truthStates(timeIdx,4);
                    true_q(count,mcIdx) = loadedFile.run(1,mcIdx).truthStates(timeIdx,5);
                    true_r(count,mcIdx) = loadedFile.run(1,mcIdx).truthStates(timeIdx,6);
                    true_lat(count,mcIdx) = loadedFile.run(1,mcIdx).truthStates(timeIdx,7);
                    true_long(count,mcIdx) = loadedFile.run(1,mcIdx).truthStates(timeIdx,8);
                    true_alt(count,mcIdx) = loadedFile.run(1,mcIdx).truthStates(timeIdx,9);
                    true_phi(count,mcIdx) = loadedFile.run(1,mcIdx).truthStates(timeIdx,10);
                    true_theta(count,mcIdx) = loadedFile.run(1,mcIdx).truthStates(timeIdx,11);
                    true_psi(count,mcIdx) = loadedFile.run(1,mcIdx).truthStates(timeIdx,12);
                    true_clkBias(count,mcIdx) = loadedFile.run(1,mcIdx).truthStates(timeIdx,13);
                    true_clkDrift(count,mcIdx) = loadedFile.run(1,mcIdx).truthStates(timeIdx,14);

                    time(count) = loadedFile.run(1,1).tout(timeIdx);
                end

                rmseU(count) = rms(u(count,mcGoodIdx) - true_u(count,mcGoodIdx),"omitnan");
                rmseV(count) = rms(v(count,mcGoodIdx) - true_v(count,mcGoodIdx),"omitnan");
                rmseW(count) = rms(w(count,mcGoodIdx) - true_w(count,mcGoodIdx),"omitnan");
                rmseSpeed(count) = rms(speed(count,mcGoodIdx) - true_speed(count,mcGoodIdx),"omitnan");
                rmseP(count) = rms(p(count,mcGoodIdx)- true_p(count,mcGoodIdx),"omitnan");
                rmseQ(count) = rms(q(count,mcGoodIdx) - true_q(count,mcGoodIdx),"omitnan");
                rmseR(count) = rms(r(count,mcGoodIdx) - true_r(count,mcGoodIdx),"omitnan");
                rmseLat(count) = rms(lat(count,mcGoodIdx) - true_lat(count,mcGoodIdx),"omitnan");
                rmseLong(count) = rms(long(count,mcGoodIdx) - true_long(count,mcGoodIdx),"omitnan");
                rmseAlt(count) = rms(alt(count,mcGoodIdx) - true_alt(count,mcGoodIdx),"omitnan");
                rmsePhi(count) = rms(phi(count,mcGoodIdx) - true_phi(count,mcGoodIdx),"omitnan");
                rmseTheta(count) = rms(theta(count,mcGoodIdx) - true_theta(count,mcGoodIdx),"omitnan");
                rmsePsi(count) = rms(psi(count,mcGoodIdx) - true_psi(count,mcGoodIdx),"omitnan");
                rmseClkBias(count) = rms(clkBias(count,mcGoodIdx) - true_clkBias(count,mcGoodIdx),"omitnan");
                rmseClkDrift(count) = rms(clkDrift(count,mcGoodIdx) - true_clkDrift(count,mcGoodIdx),"omitnan");

                stop = 1;
                count = count + 1;


            end
            obj.outputs = loadedFile.run;
            obj.timeUpdateTime = obj.outputs.tout;
            obj.measUpdateTime = obj.timeUpdateTime(1):1/50:obj.timeUpdateTime(end);
            dtMeas = diff(obj.measUpdateTime);
            dtTime = diff(obj.timeUpdateTime);

            scale = dtMeas(1)/dtTime(1);

            obj.timeIdx = 1:scale:length(obj.timeUpdateTime);

            obj.numMonteCarlos = 1;

            %% Velocity RMS
            figure('Position',[400 200 900 800])
            tiledlayout(3,1)
            nexttile
            hold on
            title('Velocity RMSE')
            fill([60 150 150 60],[0 0 30 30],[0.961 0.961 0.863],'EdgeColor','none')
            plot(time,rmseU,LineWidth=2,Color=obj.primaryColor)
            ylabel('\deltaV_N [m/s]')
            ax = gca;
            ax.FontSize = 16;
            nexttile
            hold on
            fill([60 150 150 60],[0 0 30 30],[0.961 0.961 0.863],'EdgeColor','none')
            plot(time,rmseV,LineWidth=2,Color=obj.primaryColor)
            ylabel('\deltaV_E [m/s]')
            ax = gca;
            ax.FontSize = 16;
            nexttile
            hold on
            fill([60 150 150 60],[0 0 30 30],[0.961 0.961 0.863],'EdgeColor','none')
            plot(time,rmseW,LineWidth=2,Color=obj.primaryColor)
            ylabel('\deltaV_D[m/s]')
            xlabel('Time [s]')
            ax = gca;
            ax.FontSize = 16;


            %% Angular Rates RMS
            figure('Position',[200 200 900 800])
            tiledlayout(3,1)
            nexttile
            hold on
            title('Angular Rate RMSE')
            fill([60 150 150 60],[0 0 180 180],[0.961 0.961 0.863],'EdgeColor','none')
            plot(time,rmseP*180/pi,LineWidth=2,Color=obj.primaryColor)
            ylabel('\deltap [deg/s]')
            ax = gca;
            ax.FontSize = 16;
            ylim([0 180])
            yticks([0 45 90 135 180])
            nexttile
            hold on
            fill([60 150 150 60],[0 0 180 180],[0.961 0.961 0.863],'EdgeColor','none')
            plot(time,rmseQ*180/pi,LineWidth=2,Color=obj.primaryColor)
            ylabel('\deltaq [deg/s]')
            ax = gca;
            ax.FontSize = 16;
            ylim([0 180])
            yticks([0 45 90 135 180])
            nexttile
            hold on
            fill([60 150 150 60],[0 0 180 180],[0.961 0.961 0.863],'EdgeColor','none')
            plot(time,rmseR*180/pi,LineWidth=2,Color=obj.primaryColor)
            ylabel('\deltar [deg/s]')
            xlabel('Time [s]')
            ylim([0 180])
            yticks([0 45 90 135 180])
            ax = gca;
            ax.FontSize = 16;

            %% Position RMS
            figure('Position',[600 200 900 800])
            tiledlayout(3,1)
            nexttile
            hold on
            title('Position RMSE')
            fill([60 150 150 60],[0 0 4e-3 4e-3],[0.961 0.961 0.863],'EdgeColor','none')
            plot(time,rmseLong*180/pi,LineWidth=2,Color=obj.primaryColor)
            ylabel('\delta\lambda [deg]')
            ax = gca;
            ax.FontSize = 16;
            % ylim([0 8])
            nexttile
            hold on
            fill([60 150 150 60],[0 0 4e-3 4e-3],[0.961 0.961 0.863],'EdgeColor','none')
            plot(time,rmseLat*180/pi,LineWidth=2,Color=obj.primaryColor)
            ylabel('\deltaL [deg]')
            ax = gca;
            ax.FontSize = 16;
            % ylim([0 8])
            nexttile
            hold on
            fill([60 150 150 60],[0 0 100 100],[0.961 0.961 0.863],'EdgeColor','none')
            plot(time,rmseAlt,LineWidth=2,Color=obj.primaryColor)
            ylabel('\deltah[m]')
            xlabel('Time [s]')
            ax = gca;
            ax.FontSize = 16;
            ylim([0 50])

            %% Euler Angles RMS
            figure('Position',[800 200 900 800])
            tiledlayout(3,1)
            nexttile
            hold on
            title('Euler Angles RMSE')
            fill([60 150 150 60],[0 0 180 180],[0.961 0.961 0.863],'EdgeColor','none')
            plot(time,abs(wrapTo180(rmsePhi*180/pi)),LineWidth=2,Color=obj.primaryColor)
            ylabel('\delta\phi [deg]')
            ax = gca;
            ax.FontSize = 16;
            ylim([0 180])
            yticks([0 45 90 135 180])
            nexttile
            hold on
            fill([60 150 150 60],[0 0 180 180],[0.961 0.961 0.863],'EdgeColor','none')
            plot(time,rmseTheta*180/pi,LineWidth=2,Color=obj.primaryColor)
            ylabel('\delta\theta [deg]')
            ax = gca;
            ax.FontSize = 16;
            ylim([0 180])
            yticks([0 45 90 135 180])
            nexttile
            hold on
            fill([60 150 150 60],[0 0 180 180],[0.961 0.961 0.863],'EdgeColor','none')
            plot(time,abs(wrapTo180(rmsePsi*180/pi)),LineWidth=2,Color=obj.primaryColor)
            ylabel('\delta\psi [deg]')
            xlabel('Time [s]')
            ylim([0 180])
            yticks([0 45 90 135 180])
            ax = gca;
            ax.FontSize = 16;

            %% Clock Terms
            figure('Position',[1000 200 900 800])
            tiledlayout(2,1)
            nexttile
            hold on
            title('Clock Bias/Drift RMSE')
            fill([60 150 150 60],[0 0 4 4],[0.961 0.961 0.863],'EdgeColor','none')
            plot(time,rmseClkBias,LineWidth=2,Color=obj.primaryColor)
            ylabel('Clock Bias Error [m]')
            ax = gca;
            ax.FontSize = 16;
            nexttile
            hold on
            fill([60 150 150 60],[0 0 0.15 0.15],[0.961 0.961 0.863],'EdgeColor','none')
            plot(time,rmseClkDrift,LineWidth=2,Color=obj.primaryColor)

            ylabel('Clock Drift Error [m/s]')
            ax = gca;
            ax.FontSize = 16;

            %% Code Phase and Carrier Frequency Error
            figure('Position',[0 200 900 800])
            tiledlayout(2,1)
            nexttile
            hold on
            title('Residuals RMSE')
            fill([60 150 150 60],[0 0 100 100],[0.961 0.961 0.863],'EdgeColor','none')
            plot(time,rmseCodePhase,LineWidth=2,Color=obj.primaryColor)
            ylabel('Code Phase [m]')
            ax = gca;
            ax.FontSize = 16;
            nexttile
            hold on
            fill([60 150 150 60],[0 0 3 3],[0.961 0.961 0.863],'EdgeColor','none')
            plot(time,rmseCarrFreq,LineWidth=2,Color=obj.primaryColor)
            ylabel('Carrier Frequency [m/s]')
            ax = gca;
            ax.FontSize = 16;

        end



    end
end

