%% H-Infinity Controller Design for Jiffy Jerboa in Hover Configuration
% 
%   History:
%        04.10.2021: Created, TVG
%
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

Gp_tf = zpk(CT*(s*eye(size(A))-A)^-1*B+D);

% % Extract State Space Data From Plant
% [Ag, B g, Cg, Dg] = ssdata(Gp_tf);
 
Ag = A;
Bg = B;
Cg = CT;
Dg = D;

[Ag, Bg, Cg, Dg] = minreal(Ag, Bg, Cg, Dg);

% Shift the plant to stabilize it since Hinf is undefined on jw axis
Ag = Ag - 0.05*eye(size(Ag));

%% Weights
dcgain1 = 1e-2;
mag1 = 10;
hfgain1 = 1e2; 
% multfact1 = 0.01;
multfact1 = 1;

dcgain2 = 1e-1;
mag2 = 10;
hfgain2 = 1e1;
% multfact2 = 0.1;
multfact2 = 1;

% dcgain3 = 1e-2;
% mag3 = 1e2;
% hfgain3 = 1e1;
% multfact3 = 0.001;

Wd1 = makeweight(dcgain1,mag1,hfgain1);
Wd1_tf = tf(Wd1);

Wd2 = makeweight(dcgain1,mag1,hfgain1);
Wd2_tf = tf(Wd2);

Wd3 = makeweight(dcgain1,mag1,hfgain1);
Wd3_tf = tf(Wd3);

Wd4 = makeweight(dcgain1,mag1,hfgain1);
Wd4_tf = tf(Wd4);

Hinf.Wd = multfact1*[Wd1_tf, 0, 0, 0; 0, Wd2_tf, 0, 0; 0, 0, Wd3_tf, 0; 0, 0, 0, Wd4_tf];
Wdss = ss(Hinf.Wd);

[Ad, Bd, Cd, Dd] = ssdata(Wdss);
[Ad, Bd, Cd, Dd] = minreal(Ad, Bd, Cd, Dd);

Wp1 = makeweight(dcgain2,mag2,hfgain2);
Wp1_tf = tf(Wp1);

Wp2 = makeweight(dcgain2,mag2,hfgain2);
Wp2_tf = tf(Wp2);

Wp3 = makeweight(dcgain2,mag2,hfgain2);
Wp3_tf = tf(Wp3);

Wp4 = makeweight(dcgain2,mag2,hfgain2);
Wp4_tf = tf(Wp4);

Hinf.Wp = multfact2*[Wp1_tf, 0, 0, 0; 0, Wp2_tf, 0, 0; 0, 0, Wp3_tf, 0; 0, 0, 0, Wp4_tf];

Wpss = ss(Hinf.Wp);
[Ap, Bp, Cp, Dp] = ssdata(Wpss);
[Ap, Bp, Cp, Dp] = minreal(Ap, Bp, Cp, Dp);

% Youla Weighting [8 x 8]
eps = 0.001;
Wynum = mat2cell(eps*eye(8),ones(8,1),ones(8,1)); 
Wyden = mat2cell(ones(8),ones(8,1),ones(8,1));

Hinf.Wy = tf(Wynum, Wyden);

Wyss = ss(Hinf.Wy);
[Au, Bu, Cu, Du] = ssdata(Wyss);
[Au, Bu, Cu, Du] = minreal(Au, Bu, Cu, Du);

%% Compute Augmented Plant
[A, B1, B2, C1, C2, D11, D12, D21, D22] = ...
    augss(Ag, Bg, Cg, Dg,...
          Ap, Bp, Cp, Dp,...
          Au, Bu, Cu, Du,...
          Ad, Bd, Cd, Dd);

B = [B1, B2];
C = [C1; C2];
D = [D11, D12; D21, D22];

Gaug = ss(A, B, C, D);
Hinf.Gaug = minreal(Gaug);

      
%% Compute Hinf Controller
[Kss, Cl, Gam] = hinfsyn( Hinf.Gaug, 4, 8 );
[acp, bcp, ccp, dcp] = ssdata( Kss );

acp = acp + 0.05*eye(size(acp));

Kss = ss(acp, bcp, ccp, dcp);
K = balred(tf(Kss),2);
Hinf.K = minreal(K);
fprintf("minreal(K) done \n")

Ly = Gp_tf * Hinf.K;
Hinf.Ly = minreal(Ly);
fprintf("minreal(Ly) done \n")

Ty = feedback(Hinf.Ly, eye(4));
Hinf.Ty = minreal(Ty);
fprintf("minreal(Ty) done \n")

%%
Sy = (eye(4) - Hinf.Ty);
Hinf.Sy = minreal(Sy  );
fprintf("minreal(Sy) done \n")

Y = Hinf.K * (eye(size(Hinf.Ly)) + Hinf.Ly)^(-1);
Hinf.Y = minreal(Y );
fprintf("minreal(Y) done \n")

% Hinf.Tu = minreal( Y * Gp_tf );
% fprintf("minreal(Tu) done \n")
% 
% Hinf.Su = eye(size(Hinf.Tu)) - Hinf.Tu;
% %%
% w=logspace(-3,3);
% [SVWpSy] = sigma(Hinf.Wp * Hinf.Sy, w);
% LmWpSy = 20*log10(SVWpSy(1,:));
% 
% [SVWdTy]=sigma(Hinf.Wd * Hinf.Ty,w);
% LmWdTy=20*log10(SVWdTy(1,:));
% 
% LmRp=20*log10(SVWpSy(1,:)+SVWdTy(1,:));

%% Plotting
t=0:0.01:10;
u = zeros(size(t));
u(t>0.5) = -ones();
figure('Position',[0 0 800 800])
hold on
lsim(Ty(1,1),u,t)   
lsim(Ty(2,1),u,t)  
lsim(Ty(3,1),u,t)
lsim(Ty(4,1),u,t)
lsim(Ty(5,1),u,t)
lsim(Ty(6,1),u,t)
lsim(Ty(7,1),u,t)
lsim(Ty(8,1),u,t)

