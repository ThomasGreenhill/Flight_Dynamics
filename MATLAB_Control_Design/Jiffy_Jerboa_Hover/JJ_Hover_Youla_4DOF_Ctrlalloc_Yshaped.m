%% Youla Controller Design for Jiffy Jerboa in Hover Configuration
%
%   History:
%        04.05.2021: Created, TVG
%        04.09.2021: Got responses looking good but actuator effort is high
%        04.14.2021: Tested different shapes of Ty, still no good results
%        04.22.2021: Changed to shape Youla rather than Ty
%        04.23.2021: Tuned Youla, getting good results 

%% Startup
clc;
clear;
close all;
warning('off', 'all')
mkdir './Figures'
warning('on', 'all')
addpath('./MIMO_Functions/')

%% Aircraft Inertial Properties
CG_xyz = [-3.5544829, 0, 0].';
m = 13000 / 9.81;

Ixx = 2353;
Iyy = 2327;
Izz = 3974;
Ixz = 1;

%% Trim condition is taken as zero pitch angle and zero steady-state velocities
The0 = 0;
U0 = 0;
V0 = 0;
W0 = 0;

%% Plant Model
[A, B, CT, D, E] = JJ_Hover_State_Space_4DOF(false, m, Ixx, Iyy, Izz, Ixz, The0, U0, V0, W0);

% Invert E onto A to get things in standard form
A = E \ A;
B = E \ B;

s = tf('s');
SS = ss(A, B, CT, D);
Gp_tf = CT * (s * eye(size(A)) - A)^-1 * B + D;

%% Compute the SM form of GP

TFM = cell(6, 8);
num = cell(6, 8);
den = cell(6, 8);

syms s

SI = eye(size(A)) * s;

p = simplify(CT*adjoint(SI - A)*B+D*det(SI - A));
d = 1 / simplify(det(SI - A));

Gp_sym = p .* d;
TFM_tf = minreal(zpk(sym2tf(Gp_sym)));

% Control allocation
% Leave offset vector zero for now
% c = zeros(8,1);
Gps = pinv(Gp_sym);
Gpdes = diag([1/s, 1/s, 1/s, 1/s]);
Ctrlalloc_sym = double(Gps*Gpdes);
p_alloc = p*Ctrlalloc_sym;

Gp_sym = Gp_sym*Ctrlalloc_sym;

[UL, Mp, UR] = MNsmithmcmillanForm(p_alloc, d);

%% Come up with a controller using Youla Parametrization
wn = 3.75;
zet = 0.779703;
Tz1 = 0.2;
Tp1 = 0.35;
K = wn^2;

% Shape the Youla Parameter
Y11 = K * (Tz1*s + 1) * s / ((s^2 + 2*zet*wn*s + wn^2) * (Tp1*s + 1)^2);
sigma(zpk(sym2tf(Y11)))

Y22 = Y11;
Y33 = Y11;
Y44 = Y11;

MY = sym(zeros(4));
MY(1, 1) = Y11;
MY(2, 2) = Y22;
MY(3, 3) = Y33;
MY(4, 4) = Y44;

Youla = UR * MY * UL;

Youla_tf = minreal(sym2tf(Youla));

Youlaalloc_tf = Ctrlalloc_sym*Youla_tf;

Gc_sym = simplify(Youla*(eye(size(Gp_sym * Youla)) - Gp_sym * Youla)^-1);
Gc_tf = zpk(sym2tf(Gc_sym));

Ly_tf = zpk(sym2tf(Gp_sym * Gc_sym));

Ty_tf = minreal((eye(size(Ly_tf)) + Ly_tf)^-1 * Ly_tf);
Sy_tf = minreal((eye(size(Ly_tf)) + Ly_tf)^-1);

Lu_tf = minreal(Gc_tf * Gp_tf);
I = eye(size(Lu_tf));
% Youla2 = (I + Lu_tf)^-1 * Gc_tf;
% Tu = (I + Lu_tf)^-1 * Gc_tf * Gp_tf;
% Su = (I + Gc_tf * Gp_tf)^-1;

%% Principal Gains of $L_y, S_y, T_y$
w = logspace(-1, 3, 1000);
figure('Position', [0, 0, 800, 800])
hold on
sigma(Ly_tf, w)
sigma(Sy_tf, w)
sigma(Ty_tf, w)
grid on
legend
set(findall(gcf, 'type', 'line'), 'linewidth', 2);
set(findall(gcf, 'type', 'axes'), 'fontsize', 14);
set(findall(gcf, 'type', 'legend'), 'fontsize', 14);

%% Output Step Responses to First Reference Input (Altitude Step)
t = 0:0.01:10;

figure('Position', [0, 0, 800, 800])
hold on
u = zeros(size(t));
u(t > 0.5) = -9.81*t(t > 0.5); % Hover 
lsim(Ty_tf(1, 1), u, t)
lsim(Ty_tf(2, 1), u, t)
lsim(Ty_tf(3, 1), u, t)
lsim(Ty_tf(4, 1), u, t)


grid on
legend('Response of w(t)', 'Response of p(t)', 'Response of q(t)', 'Response of r(t)')
title('Sensor Response w(t), p(t), q(t) and r(t) to Step Input at w(t))', 'FontSize', 18);

set(findall(gcf, 'type', 'line'), 'linewidth', 2);
set(findall(gcf, 'type', 'axes'), 'fontsize', 14);
set(findall(gcf, 'type', 'legend'), 'fontsize', 14);
ylim([-2,2])
saveas(gcf, "./Figures/zstep_sensor_resp.jpg")

%% Actuator Response to First Reference Input (Altitude Step)
figure('Position', [0, 0, 800, 800])
hold on
lsim(Youlaalloc_tf(1, 1), u, t)
lsim(Youlaalloc_tf(2, 1), u, t)
lsim(Youlaalloc_tf(3, 1), u, t)
lsim(Youlaalloc_tf(4, 1), u, t)
lsim(Youlaalloc_tf(5, 1), u, t)
lsim(Youlaalloc_tf(6, 1), u, t)
lsim(Youlaalloc_tf(7, 1), u, t)
lsim(Youlaalloc_tf(8, 1), u, t)

grid on
legend('Response of d_1', 'Response of d_2', 'Response of d_3', 'Response of d_4', 'Response of d_5', 'Response of d_6', 'Response of d_7', 'Response of d_8')
title('Actuator Response w(t), p(t), q(t) and r(t) to Step Input at w(t))', 'FontSize', 18);

set(findall(gcf, 'type', 'line'), 'linewidth', 2);
set(findall(gcf, 'type', 'axes'), 'fontsize', 14);
set(findall(gcf, 'type', 'legend'), 'fontsize', 14);
ylim([-2,2])
saveas(gcf, "./Figures/zstep_sensor_resp.jpg")

%% Output Step Responses to 2nd Reference Input (Roll Rate Step) from Hover
t = 0:0.01:10;

figure('Position', [0, 0, 800, 800])
hold on
u = zeros(length(t),2);
u(:,1) = -9.81*t; % Hover 
u(t > 4,2) = -(360/10)*pi/180; % 10 seconds for full yaw rotation
ind = [1,2];

lsim(Ty_tf(1, ind), u, t)
lsim(Ty_tf(2, ind), u, t)
lsim(Ty_tf(3, ind), u, t)
lsim(Ty_tf(4, ind), u, t)

grid on
legend('Response of w(t)', 'Response of p(t)', 'Response of q(t)', 'Response of r(t)')
title('Sensor Response w(t), p(t), q(t) and r(t) to Step Input at p(t))', 'FontSize', 18);

set(findall(gcf, 'type', 'line'), 'linewidth', 2);
set(findall(gcf, 'type', 'axes'), 'fontsize', 14);
set(findall(gcf, 'type', 'legend'), 'fontsize', 14);
ylim([-2,2])
saveas(gcf, "./Figures/pstep_sensor_resp.jpg")

%% Actuator Response to 2nd Reference Input (Roll Rate Step) from Hover
figure('Position', [0, 0, 800, 800])
hold on
lsim(Youlaalloc_tf(1, ind), u, t)
lsim(Youlaalloc_tf(2, ind), u, t)
lsim(Youlaalloc_tf(3, ind), u, t)
lsim(Youlaalloc_tf(4, ind), u, t)
lsim(Youlaalloc_tf(5, ind), u, t)
lsim(Youlaalloc_tf(6, ind), u, t)
lsim(Youlaalloc_tf(7, ind), u, t)
lsim(Youlaalloc_tf(8, ind), u, t)

grid on
legend('Response of d_1', 'Response of d_2', 'Response of d_3', 'Response of d_4', 'Response of d_5', 'Response of d_6', 'Response of d_7', 'Response of d_8')
title('Actuator Response w(t), p(t), q(t) and r(t) to Step Input at p(t))', 'FontSize', 18);

set(findall(gcf, 'type', 'line'), 'linewidth', 2);
set(findall(gcf, 'type', 'axes'), 'fontsize', 14);
set(findall(gcf, 'type', 'legend'), 'fontsize', 14);
ylim([-2,2])
saveas(gcf, "./Figures/pstep_actuator_resp.jpg")

%% Output Step Responses to 3rd Reference Input (Pitch Rate Step) from Hover
t = 0:0.01:10;

figure('Position', [0, 0, 800, 800])
hold on
u = zeros(length(t),2);
u(:,1) = -9.81*t; % Hover 
u(t > 4,2) = -(360/10)*pi/180; % 20 seconds for full yaw rotation
ind = [1,3];

lsim(Ty_tf(1, ind), u, t)
lsim(Ty_tf(2, ind), u, t)
lsim(Ty_tf(3, ind), u, t)
lsim(Ty_tf(4, ind), u, t)



grid on
legend('Response of w(t)', 'Response of p(t)', 'Response of q(t)', 'Response of r(t)')
title('Sensor Response w(t), p(t), q(t) and r(t) to Step Input at q(t))', 'FontSize', 18);

set(findall(gcf, 'type', 'line'), 'linewidth', 2);
set(findall(gcf, 'type', 'axes'), 'fontsize', 14);
set(findall(gcf, 'type', 'legend'), 'fontsize', 14);
ylim([-2,2])
saveas(gcf, "./Figures/qstep_sensor_resp.jpg")

%% Actuator Response to 3rd Reference Input (Pitch Rate Step) from Hover
figure('Position', [0, 0, 800, 800])
hold on
lsim(Youlaalloc_tf(1, ind), u, t)
lsim(Youlaalloc_tf(2, ind), u, t)
lsim(Youlaalloc_tf(3, ind), u, t)
lsim(Youlaalloc_tf(4, ind), u, t)
lsim(Youlaalloc_tf(5, ind), u, t)
lsim(Youlaalloc_tf(6, ind), u, t)
lsim(Youlaalloc_tf(7, ind), u, t)
lsim(Youlaalloc_tf(8, ind), u, t)

grid on
legend('Response of d_1', 'Response of d_2', 'Response of d_3', 'Response of d_4', 'Response of d_5', 'Response of d_6', 'Response of d_7', 'Response of d_8')
title('Actuator Response w(t), p(t), q(t) and r(t) to Step Input at q(t))', 'FontSize', 18);

set(findall(gcf, 'type', 'line'), 'linewidth', 2);
set(findall(gcf, 'type', 'axes'), 'fontsize', 14);
set(findall(gcf, 'type', 'legend'), 'fontsize', 14);
ylim([-2,2])
saveas(gcf, "./Figures/qstep_actuator_resp.jpg")


%% Output Step Responses to 4th Reference Input (Yaw Rate Slow Step) from Hover
t = 0:0.01:20;

figure('Position', [0, 0, 800, 800])
hold on
u = zeros(length(t),2);
u(:,1) = -9.81*t; % Hover 
ramptime = 5;
u(t > 4 & t < (4+ramptime),2) = -(360/20)*pi/180*(t(t > 4 & t < (4+ramptime))-4)/ramptime; 
u(t >= 4+ramptime,2) = -(360/20)*pi/180; % 20 seconds for full yaw rotation
ind = [1,4];

lsim(Ty_tf(1, ind), u, t)
lsim(Ty_tf(2, ind), u, t)
lsim(Ty_tf(3, ind), u, t)
lsim(Ty_tf(4, ind), u, t)



grid on
legend('Response of w(t)', 'Response of p(t)', 'Response of q(t)', 'Response of r(t)')
title('Sensor Response w(t), p(t), q(t) and r(t) to Step Input at r(t))', 'FontSize', 18);

set(findall(gcf, 'type', 'line'), 'linewidth', 2);
set(findall(gcf, 'type', 'axes'), 'fontsize', 14);
set(findall(gcf, 'type', 'legend'), 'fontsize', 14);
ylim([-2,2])
saveas(gcf, "./Figures/rstep_sensor_resp.jpg")

%% Actuator Response to 4th Reference Input (Yaw Rate Step) from Hover
figure('Position', [0, 0, 800, 800])
hold on
lsim(Youlaalloc_tf(1, ind), u, t)
lsim(Youlaalloc_tf(2, ind), u, t)
lsim(Youlaalloc_tf(3, ind), u, t)
lsim(Youlaalloc_tf(4, ind), u, t)
lsim(Youlaalloc_tf(5, ind), u, t)
lsim(Youlaalloc_tf(6, ind), u, t)
lsim(Youlaalloc_tf(7, ind), u, t)
lsim(Youlaalloc_tf(8, ind), u, t)

grid on
legend('Response of d_1', 'Response of d_2', 'Response of d_3', 'Response of d_4', 'Response of d_5', 'Response of d_6', 'Response of d_7', 'Response of d_8')
title('Actuator Response w(t), p(t), q(t) and r(t) to Step Input at r(t))', 'FontSize', 18);

set(findall(gcf, 'type', 'line'), 'linewidth', 2);
set(findall(gcf, 'type', 'axes'), 'fontsize', 14);
set(findall(gcf, 'type', 'legend'), 'fontsize', 14);
ylim([-2,2])
saveas(gcf, "./Figures/rstep_actuator_resp.jpg")


