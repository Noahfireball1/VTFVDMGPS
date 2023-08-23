function simPlot(sim,time)

%% Geoplot Trajectory
figure('units','normalized','outerposition',[0 0 1 1])
geoplot(sim.trueLAT,sim.trueLONG,'--k','LineWidth',2)
hold on
geoplot(sim.noiseLAT,sim.noiseLONG,'LineWidth',1.5)

%% States
figure('units','normalized','outerposition',[0 0 1 1])
tiledlayout(4,3)
nexttile
hold on
title('Velocity: U')
xlabel('time [s]')
ylabel('Velocity [m/s]')
plot(time,sim.truthStates(:,1),'--k','LineWidth',2)
plot(time,sim.noiseTrajectory(:,1),'LineWidth',1.5)

nexttile
hold on
title('Velocity: V')
xlabel('time [s]')
ylabel('Velocity [m/s]')
plot(time,sim.truthStates(:,2),'--k','LineWidth',2)
plot(time,sim.noiseTrajectory(:,2),'LineWidth',1.5)

nexttile
hold on
title('Velocity: W')
xlabel('time [s]')
ylabel('Velocity [m/s]')
plot(time,sim.truthStates(:,3),'--k','LineWidth',2)
plot(time,sim.noiseTrajectory(:,3),'LineWidth',1.5)

nexttile
hold on
title('Rotation Rate: phi dot')
xlabel('time [s]')
ylabel('roll rate [deg/s]')
plot(time,sim.truthStates(:,4),'--k','LineWidth',2)
plot(time,sim.noiseTrajectory(:,4),'LineWidth',1.5)

nexttile
hold on
title('Rotation Rate: theta dot')
xlabel('time [s]')
ylabel('pitch rate [deg/s]')
plot(time,sim.truthStates(:,5).*(180/pi),'--k','LineWidth',2)
plot(time,sim.noiseTrajectory(:,5).*(180/pi),'LineWidth',1.5)

nexttile
hold on
title('Rotation Rate: psi dot')
xlabel('time [s]')
ylabel('heading rate [deg/s]')
plot(time,sim.truthStates(:,6).*(180/pi),'--k','LineWidth',2)
plot(time,sim.noiseTrajectory(:,6).*(180/pi),'LineWidth',1.5)

nexttile
hold on
title('Position: X')
ylabel('Position [m]')
xlabel('time [s]')
plot(time,sim.truthStates(:,7),'--k','LineWidth',2)
plot(time,sim.noiseTrajectory(:,7),'LineWidth',1.5)

nexttile
hold on
title('Position: Y')
xlabel('time [s]')
ylabel('Position [m]')
plot(time,sim.truthStates(:,8),'--k','LineWidth',2)
plot(time,sim.noiseTrajectory(:,8),'LineWidth',1.5)

nexttile
hold on
title('Position: Z')
xlabel('time [s]')
ylabel('Position [m]')
plot(time,sim.truthStates(:,9),'--k','LineWidth',2)
plot(time,sim.noiseTrajectory(:,9),'LineWidth',1.5)

nexttile
hold on
title('Rotation: phi')
xlabel('time [s]')
ylabel('Roll Angle [deg]')
plot(time,wrapTo180(sim.truthStates(:,10).*(180/pi)),'--k','LineWidth',2)
plot(time,wrapTo180(sim.noiseTrajectory(:,10).*(180/pi)),'LineWidth',1.5)

nexttile
hold on
title('Rotation: theta')
xlabel('time [s]')
ylabel('Pitch Angle [deg]')
plot(time,wrapTo180(sim.truthStates(:,11).*(180/pi)),'--k','LineWidth',2)
plot(time,wrapTo180(sim.noiseTrajectory(:,11).*(180/pi)),'LineWidth',1.5)

nexttile
hold on
title('Rotation: psi')
xlabel('time [s]')
ylabel('Heading Angle [deg]')
plot(time,wrapTo180(sim.truthStates(:,12).*(180/pi)),'--k','LineWidth',2)
plot(time,wrapTo180(sim.noiseTrajectory(:,12).*(180/pi)),'LineWidth',1.5)


%% Orientation
figure('units','normalized','outerposition',[0 0 1 1])
tiledlayout(3,3)
nexttile
hold on
title('Groundspeed')
xlabel('time [s]')
ylabel('GS [m/s]')
plot(time,sim.truthOrientation(:,1),'--k','LineWidth',2)
plot(time,sim.noiseOrientation(:,1),'LineWidth',1.5)

nexttile
hold on
title('Vertical Velocity')
xlabel('time [s]')
ylabel('Vert Velo [m/s]')
plot(time,sim.truthOrientation(:,2),'--k','LineWidth',2)
plot(time,sim.noiseOrientation(:,2),'LineWidth',1.5)


nexttile
hold on
title('Flight Path Angle')
xlabel('time [s]')
ylabel('FPA [deg]')
plot(time,sim.truthOrientation(:,3),'--k','LineWidth',2)
plot(time,sim.noiseOrientation(:,3),'LineWidth',1.5)


nexttile
hold on
title('Heading')
xlabel('time [s]')
ylabel('HDG [deg]')
plot(time,sim.truthOrientation(:,4),'--k','LineWidth',2)
plot(time,sim.noiseOrientation(:,4),'LineWidth',1.5)


nexttile
hold on
title('Track Angle')
xlabel('time [s]')
ylabel('TRK [deg]')
plot(time,sim.truthOrientation(:,5),'--k','LineWidth',2)
plot(time,sim.noiseOrientation(:,5),'LineWidth',1.5)


nexttile
hold on
title('Crab Angle')
xlabel('time [s]')
ylabel('Crab [deg]')
plot(time,sim.truthOrientation(:,6),'--k','LineWidth',2)
plot(time,sim.noiseOrientation(:,6),'LineWidth',1.5)


nexttile
hold on
title('True Airspeed')
ylabel('TAS [m/s]')
xlabel('time [s]')
plot(time,sim.truthOrientation(:,7),'--k','LineWidth',2)
plot(time,sim.noiseOrientation(:,7),'LineWidth',1.5)


nexttile
hold on
title('Angle of Attack')
xlabel('time [s]')
ylabel('AOA [deg]')
plot(time,sim.truthOrientation(:,8),'--k','LineWidth',2)
plot(time,sim.noiseOrientation(:,8),'LineWidth',1.5)


nexttile
hold on
title('Sideslip')
xlabel('time [s]')
ylabel('Sideslip [deg]')
plot(time,sim.truthOrientation(:,9),'--k','LineWidth',2)
plot(time,sim.noiseOrientation(:,9),'LineWidth',1.5)



%% Controls

figure('units','normalized','outerposition',[0 0 1 1])
tiledlayout(3,2)
nexttile
hold on
title('Stick: Lateral')
xlabel('time [s]')
ylabel('Normalized Position')
plot(time,sim.controls(:,1),'--k','LineWidth',2)

nexttile
hold on
title('Stick: Longitudinal')
xlabel('time [s]')
ylabel('Normalized Position')
plot(time,sim.controls(:,2),'--k','LineWidth',2)

nexttile
hold on
title('Pedals')
xlabel('time [s]')
ylabel('Normalized Position')
plot(time,sim.controls(:,3),'--k','LineWidth',2)

nexttile
hold on
title('Throttle')
xlabel('time [s]')
ylabel('Normalized Position')
plot(time,sim.controls(:,4),'--k','LineWidth',2)

nexttile
hold on
title('Prop Lever')
xlabel('time [s]')
ylabel('Normalized Position')
plot(time,sim.controls(:,5),'--k','LineWidth',2)

nexttile
hold on
title('Mixture Lever')
xlabel('time [s]')
ylabel('Normalized Position')
plot(time,sim.controls(:,6),'--k','LineWidth',2)

end

