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
            for timeIdx = 1:length(loadedFile.run(1,1).tout)
                for mcIdx = 1:numMCRuns
                    try
                        u(mcIdx) = loadedFile.run(1,mcIdx).estimatedStates(timeIdx,1);
                        v(mcIdx) = loadedFile.run(1,mcIdx).estimatedStates(timeIdx,2);
                        w(mcIdx) = loadedFile.run(1,mcIdx).estimatedStates(timeIdx,3);
                        speed(mcIdx) = norm([u(mcIdx) v(mcIdx) w(mcIdx)]);
                        p(mcIdx) = loadedFile.run(1,mcIdx).estimatedStates(timeIdx,4);
                        q(mcIdx) = loadedFile.run(1,mcIdx).estimatedStates(timeIdx,5);
                        r(mcIdx) = loadedFile.run(1,mcIdx).estimatedStates(timeIdx,6);
                        lat(mcIdx) = loadedFile.run(1,mcIdx).estimatedStates(timeIdx,7);
                        long(mcIdx) = loadedFile.run(1,mcIdx).estimatedStates(timeIdx,8);
                        alt(mcIdx) = loadedFile.run(1,mcIdx).estimatedStates(timeIdx,9);
                        phi(mcIdx) = loadedFile.run(1,mcIdx).estimatedStates(timeIdx,10);
                        theta(mcIdx) = loadedFile.run(1,mcIdx).estimatedStates(timeIdx,11);
                        psi(mcIdx) = loadedFile.run(1,mcIdx).estimatedStates(timeIdx,12);
                        clkBias(mcIdx) = loadedFile.run(1,mcIdx).estimatedStates(timeIdx,13);
                        clkDrift(mcIdx) = loadedFile.run(1,mcIdx).estimatedStates(timeIdx,14);
                    catch

                        u(mcIdx) = nan;
                        v(mcIdx) = nan;
                        w(mcIdx) = nan;
                        speed(mcIdx) = nan;
                        p(mcIdx) = nan;
                        q(mcIdx) = nan;
                        r(mcIdx) = nan;
                        lat(mcIdx) = nan;
                        long(mcIdx) = nan;
                        alt(mcIdx) = nan;
                        phi(mcIdx) = nan;
                        theta(mcIdx) = nan;
                        psi(mcIdx) = nan;
                        clkBias(mcIdx) = nan;
                        clkDrift(mcIdx) = nan;
                    end

                end

                rmsU(timeIdx) = rms(u(mcIdx),"omitnan");
                rmsV(timeIdx) = rms(v(mcIdx),"omitnan");
                rmsW(timeIdx) = rms(w(mcIdx),"omitnan");
                rmsSpeed(timeIdx) = rms(speed(mcIdx),"omitnan");
                rmsP(timeIdx) = rms(p(mcIdx),"omitnan");
                rmsQ(timeIdx) = rms(q(mcIdx),"omitnan");
                rmsR(timeIdx) = rms(r(mcIdx),"omitnan");
                rmsLat(timeIdx) = rms(lat(mcIdx),"omitnan");
                rmsLong(timeIdx) = rms(long(mcIdx),"omitnan");
                rmsAlt(timeIdx) = rms(alt(mcIdx),"omitnan");
                rmsPhi(timeIdx) = rms(phi(mcIdx),"omitnan");
                rmsTheta(timeIdx) = rms(theta(mcIdx),"omitnan");
                rmsPsi(timeIdx) = rms(psi(mcIdx),"omitnan");
                rmsClkBias(timeIdx) = rms(clkBias(mcIdx),"omitnan");
                rmsClkDrift(timeIdx) = rms(clkDrift(mcIdx),"omitnan");

                stop = 1;


            end
            obj.outputs = loadedFile.run;
            obj.timeUpdateTime = obj.outputs.tout;
            obj.measUpdateTime = obj.timeUpdateTime(1):1/50:obj.timeUpdateTime(end);
            dtMeas = diff(obj.measUpdateTime);
            dtTime = diff(obj.timeUpdateTime);

            scale = dtMeas(1)/dtTime(1);

            obj.timeIdx = 1:scale:length(obj.timeUpdateTime);

            obj.numMonteCarlos = 1;

        end


    end
end

