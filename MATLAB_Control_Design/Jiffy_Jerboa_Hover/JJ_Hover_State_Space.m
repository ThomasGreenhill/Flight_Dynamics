%% State Space Model for Jiffy Jerboa in Hover Configuration
% 
%   History:
%       04.05.2021: Created and debugged, TVG
%       04.07.2021: Modified C & D for position and rate control, TVG
%
%
%% State Space Model

function [A, B, C, D, E] = JJ_Hover_State_Space(CG_xyz, m, Ixx, Iyy, Izz, Ixz, The0, U0, V0, W0)

s = @(x) sin(x);
c = @(x) cos(x);

Tmax = 2150; %N
Qmax = 155; %N.m
g = 9.81;

d1_xyz = [-0.8, 1.3, 0].';
d2_xyz = [-4.959, 5.915, -0.590].';
d3_xyz = [-4.903, 3.770, -0.590].';
d4_xyz = [-4.848, 1.620, -0.590].';
d5_xyz = [-4.848, -1.620, -0.590].';
d6_xyz = [-4.903, -3.770, -0.590].';
d7_xyz = [-4.959, -5.915, -0.590].';
d8_xyz = [-0.8, -1.3, 0].';

E = [   1       0       0       0       0       0       0       0               0        0        0       0
        0       1       0       0       0       0       0       0               0        0        0       0
        0       0       1       0       0       0       0       0               0        0        0       0
        0       0       0       1       0       0       0       0               0        0        0       0
        0       0       0       0       1       0       0       0               0        0        0       0
        0       0       0       0       0       1       0       0               0        0        0       0
        0       0       0       0       0       0       1       0               0        0        0       0
        0       0       0       0       0       0       0       1               0        0        0       -Ixz/Ixx
        0       0       0       0       0       0       0       0               1        0        0       0
        0       0       0       0       0       0       0       0               0        1        0       0
        0       0       0       0       0       0       0       0               0        0        1       0
        0       0       0       0       0       0       0       -Ixz/Izz        0        0        0       1];

A = [   0       1       0       0       0       0       0           0       0           0       0       0
        0       0       0       0       0       0       0           0       -g*c(The0)  -W0     0       V0
        0       0       0       1       0       0       0           0       0           0       0       0
        0       0       0       0       0       0       g*c(The0)   W0      0           0       0       -U0
        0       0       0       0       0       1       0           0       0           0       0       0
        0       0       0       0       0       0       0           -V0     -g*s(The0)  U0      0       0
        0       0       0       0       0       0       0           1       0           0       0       0
        0       0       0       0       0       0       0           0       0           0       0       0
        0       0       0       0       0       0       0           0       0           1       0       0
        0       0       0       0       0       0       0           0       0           0       0       0
        0       0       0       0       0       0       0           0       0           0       0       1
        0       0       0       0       0       0       0           0       0           0       0       0];
    
B = zeros(12,8);
B(6,:) = -Tmax/m*ones(1,8);
B(8,:) = [Tmax*(CG_xyz(2)-d1_xyz(2))    Tmax*(CG_xyz(2)-d2_xyz(2))  Tmax*(CG_xyz(2)-d3_xyz(2))  Tmax*(CG_xyz(2)-d4_xyz(2))  Tmax*(CG_xyz(2)-d5_xyz(2))  Tmax*(CG_xyz(2)-d6_xyz(2))  Tmax*(CG_xyz(2)-d7_xyz(2))  Tmax*(CG_xyz(2)-d8_xyz(2))]/Ixx;
B(10,:) = [Tmax*(d1_xyz(1)-CG_xyz(1))   Tmax*(d2_xyz(1)-CG_xyz(1))  Tmax*(d3_xyz(1)-CG_xyz(1))  Tmax*(d4_xyz(1)-CG_xyz(1))  Tmax*(d5_xyz(1)-CG_xyz(1))  Tmax*(d6_xyz(1)-CG_xyz(1))  Tmax*(d7_xyz(1)-CG_xyz(1))  Tmax*(d8_xyz(1)-CG_xyz(1))]/Iyy;
B(12,:) = [-Qmax    Qmax    -Qmax   Qmax    -Qmax   Qmax    -Qmax   Qmax]/Izz;

% Note that accelerations (last three rows of C and D) 
% are left out as to avoid dealing with H*xdot which affects accelerations

C = [   1   0   0   0   0   0   0           0       0           0           0       0 
        0   0   1   0   0   0   0           0       0           0           0       0
        0   0   0   0   1   0   0           0       0           0           0       0   
        0   0   0   0   0   0   0           1       0           0           0       0
        0   0   0   0   0   0   0           0       0           1           0       0
        0   0   0   0   0   0   0           0       0           0           0       1];
    
D = zeros(6,8);
    
end
