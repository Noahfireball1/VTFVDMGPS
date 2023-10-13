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

            loadedFile = load(outputFilePath);
            numMCRuns = size(loadedFile.run,2);
            % Determining RMS value for each state

            count = 1;
            for timeIdx = 1:20:length(loadedFile.run(1,1).tout)
                for mcIdx = 1:numMCRuns
                    try
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

                    catch
                        u(count,mcIdx) = nan;
                        v(count,mcIdx) = nan;
                        w(count,mcIdx) = nan;
                        speed(count,mcIdx) = nan;
                        p(count,mcIdx)= nan;
                        q(count,mcIdx) = nan;
                        r(count,mcIdx) = nan;
                        lat(count,mcIdx) = nan;
                        long(count,mcIdx) = nan;
                        alt(count,mcIdx) = nan;
                        phi(count,mcIdx) = nan;
                        theta(count,mcIdx) = nan;
                        psi(count,mcIdx) = nan;
                        clkBias(count,mcIdx) = nan;
                        clkDrift(count,mcIdx) = nan;
                        codephaseError(mcIdx,:,count) = nan(1,31,1);
                        carrFreqError(mcIdx,:,count) = nan(1,31,1);

                        true_u(count,mcIdx) = nan;
                        true_v(count,mcIdx) = nan;
                        true_w(count,mcIdx) = nan;
                        true_speed(count,mcIdx) = nan;
                        true_p(count,mcIdx) = nan;
                        true_q(count,mcIdx) = nan;
                        true_r(count,mcIdx) = nan;
                        true_lat(count,mcIdx) = nan;
                        true_long(count,mcIdx) = nan;
                        true_alt(count,mcIdx) = nan;
                        true_phi(count,mcIdx) = nan;
                        true_theta(count,mcIdx) = nan;
                        true_psi(count,mcIdx) = nan;
                        true_clkBias(count,mcIdx) = nan;
                        true_clkDrift(count,mcIdx) = nan;
                    end
                    time(count) = loadedFile.run(1,1).tout(timeIdx);
                end

                rmsU(count) = rms(u(count,:),"omitnan");
                rmsV(count) = rms(v(count,:),"omitnan");
                rmsW(count) = rms(w(count,:),"omitnan");
                rmsSpeed(count) = rms(speed(count,:),"omitnan");
                rmsP(count) = rms(p(count,:),"omitnan");
                rmsQ(count) = rms(q(count,:),"omitnan");
                rmsR(count) = rms(r(count,:),"omitnan");
                rmsLat(count) = rms(lat(count,:),"omitnan");
                rmsLong(count) = rms(long(count,:),"omitnan");
                rmsAlt(count) = rms(alt(count,:),"omitnan");
                rmsPhi(count) = rms(phi(count,:),"omitnan");
                rmsTheta(count) = rms(theta(count,:),"omitnan");
                rmsPsi(count) = rms(psi(count,:),"omitnan");
                rmsClkBias(count) = rms(clkBias(count,:),"omitnan");
                rmsClkDrift(count) = rms(clkDrift(count,:),"omitnan");

                rmseU(count) = rms(u(count,:) - true_u(count,:),"omitnan");
                rmseV(count) = rms(v(count,:) - true_v(count,:),"omitnan");
                rmseW(count) = rms(w(count,:) - true_w(count,:),"omitnan");
                rmseSpeed(count) = rms(speed(count,:) - true_speed(count,:),"omitnan");
                rmseP(count) = rms(p(count,:)- true_p(count,:),"omitnan");
                rmseQ(count) = rms(q(count,:) - true_q(count,:),"omitnan");
                rmseR(count) = rms(r(count,:) - true_r(count,:),"omitnan");
                rmseLat(count) = rms(lat(count,:) - true_lat(count,:),"omitnan");
                rmseLong(count) = rms(long(count,:) - true_long(count,:),"omitnan");
                rmseAlt(count) = rms(alt(count,:) - true_alt(count,:),"omitnan");
                rmsePhi(count) = rms(phi(count,:) - true_phi(count,:),"omitnan");
                rmseTheta(count) = rms(theta(count,:) - true_theta(count,:),"omitnan");
                rmsePsi(count) = rms(psi(count,:) - true_psi(count,:),"omitnan");
                rmseClkBias(count) = rms(clkBias(count,:) - true_clkBias(count,:),"omitnan");
                rmseClkDrift(count) = rms(clkDrift(count,:) - true_clkDrift(count,:),"omitnan");
                rmseCodePhase(count,1:31) = rms(codephaseError(:,:,count),"omitnan");
                rmseCarrFreq(count,1:31) = rms(carrFreqError(:,:,count),"omitnan");

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

