%% Check of 4DOF State-Space Function Using syms
%
%   History:
%       04.09.2021: Created, TVG
%
%% Startup
clear; close all; clc

%% Create the state space
syms CG_x CG_y CG_z m Ixx Iyy Izz Ixz The0 U0 V0 W0

CG_xyz = [CG_x, CG_y, CG_z];

[A, B, C, D, E] = JJ_Hover_State_Space_4DOF(CG_xyz, m, Ixx, Iyy, Izz, Ixz, The0, U0, V0, W0);

syms z zdot w wdot phi phidot p pdot the thedot q qdot psi psidot r rdot
xdot = [zdot, wdot, phidot, pdot, thedot, qdot, psidot, rdot].';
x = [z, w, phi, p, the, q, psi, r].';

syms d1 d2 d3 d4 d5 d6 d7 d8
u = [d1, d2, d3, d4, d5, d6, d7, d8].';

%% Turn the state space back into differential equations

eom = simplify(E*xdot == A*x + B*u)

%% T test
syms T11(s) K K11 Tz1 Tz2 zet wn Tp1 Tp2 Tp3

T11(s) = K/K11*((Tz1*s+1)*(Tz2*s+1))/((s^2+2*zet*wn*s+wn^2)*(Tp1*s+1)*(Tp2*s+1)*(Tp3*s+1)^4);

subs(T11,s,0)
simplify(subs(diff(T11,s),s,0))


%% Based on Abhi's work
syms T1(s)

T1 = K*(Tz1*s+1)*(Tz2*s+1)/((s^2+2*zet*wn*s+wn^2)*(Tp1*s+1)*(Tp2*s+1)*(Tp3*s+1)^4);
subs(T1,s,0)
simplify(subs(diff(T1,s),s,0))