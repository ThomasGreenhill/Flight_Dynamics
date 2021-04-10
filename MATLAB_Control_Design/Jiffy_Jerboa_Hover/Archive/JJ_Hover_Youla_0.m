%% Youla Controller Design for Jiffy Jerboa in Hover Configuration
% 
%   History:
%        04.05.2021: Created, TVG
%
%
%% Startup
clc; clear; close all;
addpath('/Users/thomasgreenhill/MATLAB-Drive/Toolboxes/robust');
addpath('./MIMO_Functions/');

CG_xyz = [-3.5544829, 0, 0].';
m = 13000/9.81;

Ixx = 2353;
Iyy = 2327;
Izz = 3974;
Ixz = 1;

The0 = 0;
U0 = 0;
V0 = 0;
W0 = 0;

[A, B, CT, D, E] = JJ_Hover_State_Space(CG_xyz, m, Ixx, Iyy, Izz, Ixz, The0, U0, V0, W0);

A = inv(E)*A;
B = inv(E)*B;

s = tf('s');

SS = ss(A,B,CT,D);

SI = s*eye(size(A));

Gp = CT*(SI-A)^-1*B+D;

%% Come up with a controller using Youla Parametrization

TFM = cell(6,8);
num = cell(6,8);
den = cell(6,8);

syms s

SI = eye(size(A))*s;

p = simplify(CT*adjoint(SI-A)*B + D*det(SI-A));
d = 1/simplify(det(SI - A));

% for ii = 1:6
%     for jj = 1:8
%         TFM{ii,jj} = SS(ii,jj);
%         [num{ii,jj}, den{ii,jj}] = tfdata(TFM{ii,jj});
%         num{ii,jj} = cell2mat(num{ii,jj});
%         den{ii,jj} = cell2mat(den{ii,jj});
%         
%         for kk = 1:length(num{ii,jj})
%             p(ii,jj) = p(ii,jj) + num{ii,jj}(kk)*s^(length(num{ii,jj})-kk);
%         end
%         
%         for kk = 1:length(den{ii,jj})
%             d(ii,jj) = d(ii,jj) + den{ii,jj}(kk)*s^(length(den{ii,jj})-kk);
%         end
%     end
% end

p = vpa(p,2);

TFM_sym = p./d;
TFM_tf = minreal(zpk(sym2tf(TFM_sym)));

% [UL, Mp, UR] = MNsmithmcmillanForm(A, round(B,0), CT.', D);

