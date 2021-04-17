%% State Space Model for Jiffy Jerboa in Hover Configuration
%
%   History:
%       04.05.2021: Created and debugged, TVG
%       04.07.2021: Modified C & D for position and rate control, TVG
%
%

%% State Space Model

function [A, B, C, D, E] = JJ_Hover_State_Space_4DOF(CG_xyz, m, Ixx, Iyy, Izz, Ixz, The0, U0, V0, W0)

% Assume zero steady-state roll angles
Phi0 = 0;

if isa(m, 'sym')
    symflag = true;
else
    symflag = false;
end

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

d_xyz = [d1_xyz, d2_xyz, d3_xyz, d4_xyz, d5_xyz, d6_xyz, d7_xyz, d8_xyz];

CG_xyz = zeros(3,1);

if CG_xyz == false % set the CG to be on the Center of Lift
    CG_xyz(1) = mean(d_xyz(1,:));
    CG_xyz(2) = mean(d_xyz(2,:));
    CG_xyz(3) = mean(d_xyz(3,:));
end


E = [1, 0, 0, 0, 0, 0, 0, 0; ...
    0, 1, 0, 0, 0, 0, 0, 0; ...
    0, 0, 1, 0, 0, 0, 0, 0; ...
    0, 0, 0, 1, 0, 0, 0, -Ixz / Ixx; ...
    0, 0, 0, 0, 1, 0, 0, 0; ...
    0, 0, 0, 0, 0, 1, 0, 0; ...
    0, 0, 0, 0, 0, 0, 1, 0; ...
    0, 0, 0, -Ixz / Izz, 0, 0, 0, 1];

A = [0, 1, 0, 0, 0, 0, 0, 0; ...
    0, 0, -g * s(Phi0) * c(The0), -V0, g * s(The0) * c(Phi0), U0, 0, 0; ...
    0, 0, 0, 1, 0, 0, 0, 0; ...
    0, 0, 0, 0, 0, 0, 0, 0; ...
    0, 0, 0, 0, 0, 1, 0, 0; ...
    0, 0, 0, 0, 0, 0, 0, 0; ...
    0, 0, 0, 0, 0, 0, 0, 1; ...
    0, 0, 0, 0, 0, 0, 0, 0];

if symflag
    B = sym(zeros(8, 8));
else
    B = zeros(8, 8);
end

B(2, :) = -Tmax / m * ones(1, 8);
B(4, :) = Tmax / Ixx * [(CG_xyz(2) - d1_xyz(2)), (CG_xyz(2) - d2_xyz(2)), (CG_xyz(2) - d3_xyz(2)), (CG_xyz(2) - d4_xyz(2)), (CG_xyz(2) - d5_xyz(2)), (CG_xyz(2) - d6_xyz(2)), (CG_xyz(2) - d7_xyz(2)), (CG_xyz(2) - d8_xyz(2))];
B(6, :) = Tmax / Iyy * [(d1_xyz(1) - CG_xyz(1)), (d2_xyz(1) - CG_xyz(1)), (d3_xyz(1) - CG_xyz(1)), (d4_xyz(1) - CG_xyz(1)), (d5_xyz(1) - CG_xyz(1)), (d6_xyz(1) - CG_xyz(1)), (d7_xyz(1) - CG_xyz(1)), (d8_xyz(1) - CG_xyz(1))];
B(8, :) = 1 / Izz * [-Qmax, Qmax, -Qmax, Qmax, -Qmax, Qmax, -Qmax, Qmax];

% Note that accelerations (last three rows of C and D)
% are left out as to avoid dealing with H*xdot which affects accelerations

C = [1, 0, 0, 0, 0, 0, 0, 0; ...
    0, 0, 0, 1, 0, 0, 0, 0; ...
    0, 0, 0, 0, 0, 1, 0, 0; ...
    0, 0, 0, 0, 0, 0, 0, 1];

D = zeros(4, 8);

end
