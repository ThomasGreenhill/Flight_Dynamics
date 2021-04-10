%% Youla Controller Design for Jiffy Jerboa in Hover Configuration
% 
%   History:
%        04.05.2021: Created, TVG
%
%
%% Startup
clc; clear; close all;
warning('off', 'all')
mkdir './Figures'
warning('on', 'all')

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
[A, B, CT, D, E] = JJ_Hover_State_Space(CG_xyz, m, Ixx, Iyy, Izz, Ixz, The0, U0, V0, W0);

% Invert E onto A to get things in standard form
A = inv(E)*A;
B = inv(E)*B;

s = tf('s');
SS = ss(A,B,CT,D);
Gp_tf = CT*(s*eye(size(A))-A)^-1*B+D;

%% Come up with a controller using Youla Parametrization

TFM = cell(6,8);
num = cell(6,8);
den = cell(6,8);

syms s

SI = eye(size(A))*s;

p = simplify(CT*adjoint(SI-A)*B + D*det(SI-A));
d = 1/simplify(det(SI - A));

p = vpa(p,2);

TFM_sym = p./d;
TFM_tf = minreal(zpk(sym2tf(TFM_sym)));

[UL, Mp, UR] = MNsmithmcmillanForm(A, round(B,0), CT.', D);

