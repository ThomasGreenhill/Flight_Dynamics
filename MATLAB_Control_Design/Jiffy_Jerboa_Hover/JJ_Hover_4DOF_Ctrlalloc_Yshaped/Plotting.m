%% Plotting for Dynamics_Control of JJ_Hover_4DOF_Ctrlalloc_Yshaped
%
%   History:
%       04.23.2021: Creatd based on JJ_Hover_Youla_4DOF_Ctrlalloc_Yshaped.m
%
%   Notes:
%       1. Youla is shaped, not Ty
%       2. Control is directly allocated
%       3. Dynamics are based on state space, so everything is LTI
%

%% Startup
clc;
clear;
close all;
warning('off', 'all')
mkdir('./Figures')
warning('on', 'all')
addpath('../MIMO_Functions/')
formatlatex

%% Aircraft Inertial Properties
g = 9.81;
Tmax = 2150; % Max thrust per rotor (N) 
% Set CG_xyz on center of lift for now:
CG = false;
% CG_xyz = [-3.5544829, 0, 0].';
Inertia.m = 13000 / 9.81;

Inertia.Ixx = 2353;
Inertia.Iyy = 2327;
Inertia.Izz = 3974;
Inertia.Ixz = 1;

%% Trim condition is taken as zero pitch angle and zero steady-state velocities
Trim.The0 = 0;
Trim.U0 = 0;
Trim.V0 = 0;
Trim.W0 = 0;

%% Create the dynamics and controller
[Gp_sym, Gc_sym, Ctrlalloc, Youlaalloc, Ly_tf, Lu_tf, Sy_tf, Su_tf, Ty_tf, Tu_tf] = Dynamics_Control(Inertia, CG, Trim);

%% Principal Gains of L_y, S_y, T_y
w = logspace(-2, 2, 1000);
figure('Position', [0, 0, 800, 600])
hold on
sigma(Ly_tf, w)
sigma(Sy_tf, w)
sigma(Ty_tf, w)
grid on
legend
set(findall(gcf, 'type', 'line'), 'linewidth', 2);
set(findall(gcf, 'type', 'axes'), 'fontsize', 14);
set(findall(gcf, 'type', 'legend'), 'fontsize', 14);

%% Output Step Responses to First Reference Input (Altitude Position Step While Controlling Velocity)
t = 0:0.01:20;
u = zeros(size(t));
baseline = (Inertia.m*g/8)/Tmax; 
du = u;
itime = 0; % Initiate the maneuever itime seconds after t0
u(t > itime) = -g*t(t > itime); % negative is upwards, so this results in Hover 
du(t > itime) = -g; 
state_resp = lsim(Ty_tf(:, 1), u, t);
actuator_resp = lsim(Youlaalloc(:, 1), u, t);

figure('Position', [0, 0, 1000, 1000])
sgtitle(["System Time-Response to Altitude Hold Input; Craft Initially ","Supported, Support Instantaneously Removed at $t=0$"],'fontsize',30)
subplot(2,1,1)
hold on
plot(t,state_resp./t.') % System Response
plot(t,du,'k-.') % Reference Signal
title("Response of Complimentary Sensitivity Transfer Function ($T_y$)")
xlabel("Time (s)")
ylabel("Response of $\dot{\mathbf{y}}$ (Output Derivative)")
legend("$\dot w$ (z-acceleration) (m.s${^{-2}}$)", "$\dot p$ (rad.s${^{-2}}$)", "$\dot q$ (rad.s${^{-2}}$)", "$\dot r$ (rad.s${^{-2}}$)","Reference Signal $u$",'Location','best')

subplot(2,1,2)
hold on
plot(t,actuator_resp) % System Response
plot(t,baseline*ones(size(t)),'k-.') % Baseline for Hover from force balance 
title("Response of Youla Parameter ($Y$)")
xlabel("Time (s)")
ylabel("Actuator Response ($\%$ Throttle)")
legend("$\delta_1$", "$\delta_2$", "$\delta_3$", "$\delta_4$", "$\delta_5$", "$\delta_6$", "$\delta_7$", "$\delta_8$",strcat("Baseline for Hover",string(newline),"from $\Sigma F_z = 0$"),'Location','best')
saveas(gcf,"./Figures/wdot_resp.jpg")

%% Output Step Responses to Second Reference Input (Roll Angular Rate Step from Hover Initially)
ind = [1,2]; % Controls are applied to 1st and 2nd reference inputs
t = 0:0.01:10;
baseline = (Inertia.m*g/8)/Tmax; 
u = zeros(length(t),2);
du = zeros(length(t));
itime = 5; % Initiate the maneuever itime seconds after t0
u(:,1) = -g*t; % Hover 
u(t > itime,2) = -(360/10)*pi/180; % 10 seconds for full roll rotation 
du = -g; 
state_resp = lsim(Ty_tf(:, ind), u, t);
actuator_resp = lsim(Youlaalloc(:, ind), u, t);

figure('Position', [0, 0, 1000, 1000])
sgtitle(["System Time-Response to Roll Rate Input; Craft Initially","Supported, Support Instantaneously Removed at $t=0$"],'fontsize',30)
subplot(2,1,1)
hold on
plot(t,state_resp(:,1)./t.') % System Response
plot(t,state_resp(:,2:end))
plot(t,du,'k-.') % Reference Signal
plot(t,u(:,2),'k-.')
title("Response of Complimentary Sensitivity Transfer Function ($T_y$)")
xlabel("Time (s)")
ylabel("Response of $\mathbf{y}$ (Output)")
legend("$\dot w$ (z-acceleration) (m.s${^{-2}}$)", "$p$ (rad.s${^{-1}}$)", "$q$ (rad.s${^{-1}}$)", "$r$ (rad.s${^{-1}}$)","Reference Signal $u$",'Location','best')

subplot(2,1,2)
hold on
plot(t,actuator_resp) % System Response
plot(t,baseline*ones(size(t)),'k-.') % Baseline for Hover from force balance 
title("Response of Youla Parameter ($Y$)")
xlabel("Time (s)")
ylabel("Actuator Response ($\%$ Throttle)")
legend("$\delta_1$", "$\delta_2$", "$\delta_3$", "$\delta_4$", "$\delta_5$", "$\delta_6$", "$\delta_7$", "$\delta_8$",strcat("Baseline for Hover",string(newline),"from $\Sigma F_z = 0$"),'Location','best')

saveas(gcf,"./Figures/p_resp.jpg")

%% Output Step Responses to Third Reference Input (Pitch Angular Rate Step from Hover Initially)
ind = [1,3]; % Controls are applied to 1st and 2nd reference inputs
t = 0:0.01:10;
baseline = (Inertia.m*g/8)/Tmax; 
u = zeros(length(t),2);
du = zeros(length(t));
itime = 5; % Initiate the maneuever itime seconds after t0
u(:,1) = -g*t; % Hover 
u(t > itime,2) = -(360/10)*pi/180; % 10 seconds for full pitch rotation
du = -g; 
state_resp = lsim(Ty_tf(:, ind), u, t);
actuator_resp = lsim(Youlaalloc(:, ind), u, t);

figure('Position', [0, 0, 1000, 1000])
sgtitle(["System Time-Response to Pitch Rate Input; Craft Initially ","Supported, Support Instantaneously Removed at $t=0$"],'fontsize',30)
subplot(2,1,1)
hold on
plot(t,state_resp(:,1)./t.') % System Response
plot(t,state_resp(:,2:end))
plot(t,du,'k-.') % Reference Signal
plot(t,u(:,2),'k-.')
title("Response of Complimentary Sensitivity Transfer Function ($T_y$)")
xlabel("Time (s)")
ylabel("Response of $\mathbf{y}$ (Output)")
legend("$\dot w$ (z-acceleration) (m.s${^{-2}}$)", "$p$ (rad.s${^{-1}}$)", "$q$ (rad.s${^{-1}}$)", "$r$ (rad.s${^{-1}}$)","Reference Signal $u$",'Location','best')

subplot(2,1,2)
hold on
plot(t,actuator_resp) % System Response
plot(t,baseline*ones(size(t)),'k-.') % Baseline for Hover from force balance 
title("Response of Youla Parameter ($Y$)")
xlabel("Time (s)")
ylabel("Actuator Response ($\%$ Throttle)")
legend("$\delta_1$", "$\delta_2$", "$\delta_3$", "$\delta_4$", "$\delta_5$", "$\delta_6$", "$\delta_7$", "$\delta_8$",strcat("Baseline for Hover",string(newline),"from $\Sigma F_z = 0$"),'Location','best')

saveas(gcf,"./Figures/q_resp.jpg")

%% Output Step Responses to Fourth Reference Input (Yaw Angular Rate Step from Hover Initially)
ind = [1,4]; % Controls are applied to 1st and 2nd reference inputs
t = 0:0.01:20;
baseline = (Inertia.m*g/8)/Tmax; 
u = zeros(length(t),2);
du = zeros(length(t));
itime = 5; % Initiate the maneuever itime seconds after t0
u(:,1) = -g*t; % Hover 
ramptime = 5;
u(t > 4 & t < (4+ramptime),2) = -(360/20)*pi/180*(t(t > 4 & t < (4+ramptime))-4)/ramptime; 
u(t >= 4+ramptime,2) = -(360/20)*pi/180; % 20 seconds for full yaw rotation
du = -g; 
state_resp = lsim(Ty_tf(:, ind), u, t);
actuator_resp = lsim(Youlaalloc(:, ind), u, t);

figure('Position', [0, 0, 1000, 1000])
sgtitle(["System Time-Response to Yaw Rate Input; Craft Initially ","Supported, Support Instantaneously Removed at $t=0$"],'fontsize',30)
subplot(2,1,1)
hold on
plot(t,state_resp(:,1)./t.') % System Response
plot(t,state_resp(:,2:end))
plot(t,du,'k-.') % Reference Signal
plot(t,u(:,2),'k-.')
title("Response of Complimentary Sensitivity Transfer Function ($T_y$)")
xlabel("Time (s)")
ylabel("Response of $\mathbf{y}$ (Output)")
legend("$\dot w$ (z-acceleration) (m.s${^{-2}}$)", "$p$ (rad.s${^{-1}}$)", "$q$ (rad.s${^{-1}}$)", "$r$ (rad.s${^{-1}}$)","Reference Signal $u$",'Location','best')

subplot(2,1,2)
hold on
plot(t,actuator_resp) % System Response
plot(t,baseline*ones(size(t)),'k-.') % Baseline for Hover from force balance 
title("Response of Youla Parameter ($Y$)")
xlabel("Time (s)")
ylabel("Actuator Response ($\%$ Throttle)")
legend("$\delta_1$", "$\delta_2$", "$\delta_3$", "$\delta_4$", "$\delta_5$", "$\delta_6$", "$\delta_7$", "$\delta_8$",strcat("Baseline for Hover",string(newline),"from $\Sigma F_z = 0$"),'Location','best')

saveas(gcf,"./Figures/r_resp.jpg")

%% formatlatex FUNCTION DEFINITION
function formatlatex

reset(groot)
set(groot, 'defaulttextinterpreter', 'latex')
set(groot, 'defaultcolorbarticklabelinterpreter', 'latex')
set(groot, 'defaultfigureposition', [100, 100, 1000, 1200])
set(groot, 'defaultaxesticklabelinterpreter', 'latex')
set(groot, 'defaultlegendinterpreter', 'latex')
set(groot, 'defaultaxesfontsize', 18)
set(groot, 'defaultaxeslinewidth', 1)
set(groot, 'defaultscattermarker', 'o')
set(groot, 'defaultscatterlinewidth', 2)
set(groot, 'defaultlinemarkersize', 15)
set(groot, 'defaultlinelinewidth', 2.5)
set(groot, 'defaultAxesXgrid', 'on', 'defaultAxesYgrid', 'on', 'defaultAxesZgrid', 'on')
set(groot, 'defaultAxesGridLineStyle', '-.')
set(groot, 'defaultAxesXlim', [0, 2 * pi]);
set(groot, 'defaultAxesYlim', [-0.6, 0.6]);
end