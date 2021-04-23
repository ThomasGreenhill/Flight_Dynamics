%% Youla Controller Design Function for Jiffy Jerboa in Hover Configuration
%
%   History:
%        04.23.2021: Converted routine from
%           JJ_Hover_Youla_4DOF_Ctrlalloc_Yshaped.m to this function
%
%   Inputs:
%
%
%   Outputs:
%
%
%   Notes:

function [Gp_sym, Gc_sym, Ctrlalloc, Youlaalloc_tf, Ly_tf, Lualloc_tf, Sy_tf, Sualloc_tf, Ty_tf, Tualloc_tf] = Dynamics_Control(Inertia, CG, Trim)

%% Unpack inputs
m = Inertia.m;
Ixx = Inertia.Ixx;
Iyy = Inertia.Iyy;
Izz = Inertia.Izz;
Ixz = Inertia.Ixz;

The0 = Trim.The0;
U0 = Trim.U0;
V0 = Trim.V0;
W0 = Trim.W0;

%% Plant Model
[A, B, CT, D, E] = JJ_Hover_State_Space_4DOF(CG, m, Ixx, Iyy, Izz, Ixz, The0, U0, V0, W0);

% Invert E onto A to get things in standard form
A = E \ A;
B = E \ B;

s = tf('s');
Gp_tf = CT * (s * eye(size(A)) - A)^-1 * B + D;

%% Compute the SM form of GP
syms s

SI = eye(size(A)) * s;

p = simplify(CT*adjoint(SI - A)*B+D*det(SI - A));
d = 1 / simplify(det(SI - A));

Gp_sym = p .* d;

% Control allocation
% Leave offset vector zero for now
% c = zeros(8,1);
Gps = pinv(Gp_sym);
Gpdes = diag([1 / s, 1 / s, 1 / s, 1 / s]);
Ctrlalloc = double(Gps*Gpdes);
p_alloc = p * Ctrlalloc;

% Slap the control allocator on the plant
Gpalloc_sym = Gp_sym * Ctrlalloc;
Gpalloc_tf = Gp_tf * Ctrlalloc;

[UL, Mp, UR] = MNsmithmcmillanForm(p_alloc, d);

%% Come up with a controller using Youla Parametrization
wn = 3.75;
zet = 0.779703;
Tz1 = 0.2;
Tp1 = 0.35;
K = wn^2;

% Shape the Youla Parameter
Y11 = K * (Tz1 * s + 1) * s / ((s^2 + 2 * zet * wn * s + wn^2) * (Tp1 * s + 1)^2);

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

Youlaalloc_tf = Ctrlalloc*Youla_tf;

Gc_sym = simplify(Youla*(eye(size(Gpalloc_sym * Youla)) - Gpalloc_sym * Youla)^-1);
Gc_tf = zpk(sym2tf(Gc_sym));

Ly_tf = zpk(sym2tf(Gpalloc_sym * Gc_sym));

Ty_tf = minreal((eye(size(Ly_tf))+Ly_tf)^-1*Ly_tf);
Sy_tf = minreal((eye(size(Ly_tf))+Ly_tf)^-1);

Lualloc_tf = minreal(Gc_tf*Gpalloc_tf);
I = eye(size(Lualloc_tf));
Tualloc_tf = (I + Lualloc_tf)^-1 * Gc_tf * Gpalloc_tf;
Sualloc_tf = (I + Gc_tf * Gpalloc_tf)^-1;

end
