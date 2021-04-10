%% Youla Controller Design for Jiffy Jerboa in Hover Configuration
% 
%   History:
%        04.05.2021: Created, TVG
%        04.09.2021: Got responses looking good but actuator effort is high
%
%% Startup
clc; clear; close all;
warning('off', 'all')
mkdir './Figures'
warning('on', 'all')
addpath('./MIMO_Functions/')

%% Aircraft Inertial Properties
CG_xyz = [-3.5544829, 0, 0].';
m = 13000/9.81;

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
[A, B, CT, D, E] = JJ_Hover_State_Space_4DOF(CG_xyz, m, Ixx, Iyy, Izz, Ixz, The0, U0, V0, W0);

% Invert E onto A to get things in standard form
A = E\A;
B = E\B;

s = tf('s');
SS = ss(A,B,CT,D);
Gp_tf = CT*(s*eye(size(A))-A)^-1*B+D;

%% Compute the SM form of GP

TFM = cell(6,8);
num = cell(6,8);
den = cell(6,8);

syms s

SI = eye(size(A))*s;

p = simplify(CT*adjoint(SI-A)*B + D*det(SI-A));
d = 1/simplify(det(SI - A));

Gp_sym = p.*d;
vpa(Gp_sym,2)
TFM_tf = minreal(zpk(sym2tf(Gp_sym)));

[UL, Mp, UR] = MNsmithmcmillanForm(A, B, CT.', D);

%% Come up with a controller using Youla Parametrization

K11 = Mp(1,1)*s^2;

wn = 3.75; 
K = wn^2*K11;
zet = sqrt(2);
Tz1 = 2;
Tz2 = 3.5;
Tp2 = 0.5;
Tp3 = 0.5;

syms T11(s) Tp1
T11(s) = K/K11*((Tz1*s+1)*(Tz2*s+1))/((s^2+2*zet*wn*s+wn^2)*(Tp1*s+1)*(Tp2*s+1));

% Check interpolation conditions:
subs(T11,s,0)
Tp1 = double(solve(simplify(subs(diff(T11,s),s,0))==0,Tp1))

T11 = K/K11*((Tz1*s+1)*(Tz2*s+1))/((s^2+2*zet*wn*s+wn^2)*(Tp1*s+1)*(Tp2*s+1)*(Tp3*s+1)^2);
T22 = T11;
T33 = T11;
T44 = T11; 

MY = sym(zeros(8,4));
MY(1,1) = T11 / Mp(1,1);
MY(2,2) = T22 / Mp(2,2);
MY(3,3) = T33 / Mp(3,3);
MY(4,4) = T44 / Mp(4,4);

Youla = UR*MY*UL;

Youla_tf = sym2tf(Youla);

Gc_sym = simplify(Youla*(eye(size(Gp_sym*Youla))-Gp_sym*Youla)^-1);
Gc_tf = sym2tf(Gc_sym);

Ly_tf = sym2tf(Gp_sym*Gc_sym);

Ty_tf = (eye(size(Ly_tf)) + Ly_tf)^-1*Ly_tf;
Sy_tf = (eye(size(Ly_tf)) + Ly_tf)^-1;

Lu_tf = Gc_tf*Gp_tf;
I = eye(size(Lu_tf));
Youla2 = (I + Lu_tf)^-1*Gc_tf;
Tu = (I+Lu_tf)^-1*Gc_tf*Gp_tf;
Su = (I+Gc_tf*Gp_tf)^-1;

% Principal Gains of $L_y, S_y, T_y$
w = logspace(-1,3,1000);
figure('Position',[0 0 800 800])
hold on
sigma(Ly_tf, w)
sigma(Sy_tf, w)
sigma(Ty_tf, w)
grid on
legend
set(findall(gcf,'type','line'),'linewidth',2);
set(findall(gcf,'type','axes'),'fontsize',14);
set(findall(gcf,'type','legend'),'fontsize',14);

%% Output Step Responses to First Reference Input (Altitude Step)
t=0:0.01:10;

figure('Position',[0 0 800 800])
hold on
u = zeros(size(t));
u(t>0.5) = -ones();
lsim(Ty_tf(1,1),u,t)   
lsim(Ty_tf(2,1),u,t)  
lsim(Ty_tf(3,1),u,t)
lsim(Ty_tf(4,1),u,t)


grid on
legend('Response of u_1','Response of u_2','Response of u_3','Response of u_4','Response of u_5','Response of u_6','Response of u_7','Response of u_8')
title('Sensor Response z(t), p(t), q(t) and r(t) to Step Input at z(t))','FontSize',18);

set(findall(gcf,'type','line'),'linewidth',2);
set(findall(gcf,'type','axes'),'fontsize',14);
set(findall(gcf,'type','legend'),'fontsize',14);
saveas(gcf,"./Figures/z_resp.jpg")

%% Actuator Response to First Reference Input (Altitude Step)
figure('Position',[0 0 800 800])
hold on
u = zeros(size(t));
u(t>0.5) = -ones();
lsim(Youla_tf(1,1),u,t)   
lsim(Youla_tf(2,1),u,t)  
lsim(Youla_tf(3,1),u,t)
lsim(Youla_tf(4,1),u,t)
lsim(Youla_tf(5,1),u,t)
lsim(Youla_tf(6,1),u,t)
lsim(Youla_tf(7,1),u,t)
lsim(Youla_tf(8,1),u,t)

grid on
legend('Response of u_1','Response of u_2','Response of u_3','Response of u_4','Response of u_5','Response of u_6','Response of u_7','Response of u_8')
title('Sensor Response z(t), p(t), q(t) and r(t) to Step Input at z(t))','FontSize',18);

set(findall(gcf,'type','line'),'linewidth',2);
set(findall(gcf,'type','axes'),'fontsize',14);
set(findall(gcf,'type','legend'),'fontsize',14);
saveas(gcf,"./Figures/p_resp.jpg")
