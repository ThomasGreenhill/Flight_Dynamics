function Hinf = Produce_Lateral_Hinf_2( Gp )
%% Aircraft State-Space Representation Conversion
% Thomas Greenhill
% Blake Christierson
% June 8, 2020
%
% MAE 273B SQ2020
% Professor Francis Assadian

%% Notes:
%{ 

Modified and formatted from "InclassExample.m" provided by Professor 
    Assadian on June 4, 2020 

%}

%% Main Function
% Setup Variables
w = logspace(-2,2);
s = tf('s');

% Extract State Space Data From Plant
[Ag, Bg, Cg, Dg] = ssdata(Gp);
 
[Ag, Bg, Cg, Dg] = minreal(Ag, Bg, Cg, Dg);

Ag = Ag - 0.05*eye(size(Ag));

%% Weight Construction
% Complimentary Sensitivity Weighting [3 x 3]
Wd1 = makeweight(0.5,10,50); 
% Wd1 = makeweight(0.99,10,50); % Attempt for Full Plant, Not Rounded
Wd1tf = tf(Wd1);

Wd2 = makeweight(0.5,10,50); 
%Wd2 = makeweight(0.99,10,50); % Attempt for Full Plant, Not Rounded 
Wd2tf = tf(Wd2);

Wd3 = makeweight(0.5,10,50);
%Wd3 = makeweight(0.99,10,50);  % Attempt for Full Plant, Not Rounded
Wd3tf = tf(Wd3);

Hinf.Wd = [Wd1tf 0 0; 0 Wd2tf 0; 0 0 Wd3tf];

Wdss = ss(Hinf.Wd);
[Ad, Bd, Cd, Dd] = ssdata(Wdss);
[Ad, Bd, Cd, Dd] = minreal(Ad, Bd, Cd, Dd);

% Sensitivity Weighting [3 x 3]
Wp1 = makeweight(100,0.8,0.5);
Wp1tf = tf(Wp1);

Wp2 = makeweight(100,0.8,0.5);
Wp2tf = tf(Wp2);

Wp3 = makeweight(100,0.8,0.5);
Wp3tf = tf(Wp3);

Hinf.Wp = [Wp1tf 0 0;0 Wp2tf 0; 0 0 Wp3tf];

Wpss = ss(Hinf.Wp);
[Ap, Bp, Cp, Dp] = ssdata(Wpss);
[Ap, Bp, Cp, Dp] = minreal(Ap, Bp, Cp, Dp);

% Youla Weighting [2 x 2]
eps = 0.001;
Wynum = {eps 0; 0 eps}; 
Wyden = {1 1; 1 1};

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
% [acp, bcp, ccp, dcp, acl, bcl, ccl, dcl] = ...
%      hinf(A, B1, B2, C1, C2, D11, D12, D21, D22);

%[Hinf.Gam, acp, bcp, ccp, dcp, acl, bcl, ccl, dcl] = ...
%    hinfopt(A, B1, B2, C1, C2, D11, D12, D21, D22);

[Kss, Cl, Gam] = hinfsyn( Hinf.Gaug, 3, 2 );
[acp, bcp, ccp, dcp] = ssdata( Kss );

acp = acp + 0.05*eye(size(acp));

Kss = ss(acp, bcp, ccp, dcp);
K = tf(Kss);
Hinf.K = minreal(K);

Ly = Gp * Hinf.K;
Hinf.Ly = minreal(Ly);

Ty = feedback(Hinf.Ly, eye(3));
Hinf.Ty = minreal(Ty);

Sy = (eye(3) - Hinf.Ty);
Hinf.Sy = minreal(Sy  );

Y = Hinf.K * (eye(size(Hinf.Ly)) + Hinf.Ly)^(-1);
Hinf.Y = minreal(Y );

Hinf.Tu = minreal( Y * Gp );

Hinf.Su = eye(size(Hinf.Tu)) - Hinf.Tu;

% PrincipleGainAnalysis( Gp, Hinf.K, Hinf.Ly, ...
%    Hinf.Ty, Hinf.Sy, Hinf.Su, Hinf.Y, logspace(-2, 2) )

[SVWpSy] = sigma(Hinf.Wp * Hinf.Sy, w);
LmWpSy = 20*log10(SVWpSy(1,:));

[SVWdTy]=sigma(Hinf.Wd * Hinf.Ty,w);
LmWdTy=20*log10(SVWdTy(1,:));

LmRp=20*log10(SVWpSy(1,:)+SVWdTy(1,:));

figure(1)
sigma(Hinf.Ly, w);
hold on
sigma(Hinf.Wp, w);
sigma(1 / Hinf.Wd, w);
title('Open Loop and Filter Transfer Function Singular Values');
h = legend({'$L_y$','$W_p$','$1/W_{\Delta}$'}, 'interpreter', 'latex');
set(h,'string',{'$L_y$','$W_p$','$1/W_{\Delta}$'});
grid
set(gca,'Color','w')


figure(2)
sigma(Hinf.Sy,w);
hold on
sigma(Hinf.Ty,w);
hold on
sigma(Hinf.Wp,w);
hold on
sigma(1/ Hinf.Wd,w);
title('Closed Loop and Filter Transfer Function Singular Values');
h = legend({'$S_y$','$T_y$','$W_p$','$1/W_{\Delta}$'}, 'interpreter', 'latex');
set(h,'string',{'$S_y$','$T_y$','$W_p$','$1/W_{\Delta}$'});
grid
set(gca,'Color','w')


figure(3)
sigma(Hinf.Y,w)
hold on
sigma(Gp,w)
hold on
sigma(Hinf.K,w)
title('Actuator, Plant, and Controller Singular Values','interpreter', 'tex');
h = legend({'$Y$','$G_p$','$K$'}, 'interpreter', 'latex');
set(h,'string',{'$Y$','$G_p$','$K$'});
grid
set(gca,'Color','w')


figure(4)
subplot(2,1,1)
[y,t] = step(Hinf.Ty);
k=1:length(t);
plot(t,y(k,1,1),'-');
title('First Input Step Output Responses','interpreter', 'tex')
hold on
plot(t,y(k,2,1),'--')
plot(t,y(k,3,1),'-.')
legend('$\phi(t)$','$\psi(t)$','$\beta(t)$');
xlabel('Time (s)','interpreter', 'tex')
grid
set(gca,'Color','w')

subplot(2,1,2)
[y,t] = step(Hinf.Ty);
k=1:length(t);
plot(t,y(k,1,2),'-');
title('Second Input Step Output Responses','interpreter', 'tex')
hold on
plot(t,y(k,2,2),'--')
plot(t,y(k,3,2),'-.')
legend('$\phi(t)$','$\psi(t)$','$\beta(t)$');
xlabel('Time (s)','interpreter', 'tex')
grid
set(gca,'Color','w')


figure(5)
subplot(2,1,1)
t=0:0.001:0.2;
[ya] =step(Y,t);
k=1:length(t);
plot(t,ya(k,1,1),'-');
title('First Input Step Actuator Responses','interpreter', 'tex')
hold on
plot(t,ya(k,2,1),'-')
h = legend('$\delta_{a}(t)$','$\delta_{r}(t)$');
set(h,'string',{'$\delta_{a}(t)$','$\delta_{r}(t)$'});
xlabel('Time (s)','interpreter', 'tex')
grid
set(gca,'Color','w')

subplot(2,1,2)
t=0:0.001:0.2;
[ya] =step(Y,t);
k=1:length(t);
plot(t,ya(k,1,2),'-');
title('Second Input Step Actuator Responses','interpreter', 'tex')
hold on
plot(t,ya(k,2,2),'--')
h = legend('$\delta_{a}(t)$','$\delta_{r}(t)$');
set(h,'string',{'$\delta_{a}(t)$','$\delta_{r}(t)$'});
xlabel('Time (s)','interpreter', 'tex')
grid
set(gca,'Color','w')

figure(6);
semilogx(w,LmWpSy,'r',w,LmWdTy,'g',w,LmRp,'b')
grid;
title('Nominal Performance, Robust Stability, and Robust Performance','interpreter', 'tex')
h = legend('$\sigma(W_{p}*S_{y})$','$\sigma(W_{\Delta}*T_{y})$','$\sigma(W_{p}*S_{y})+\sigma(W_{\Delta}*T_{y})$');
set(h,'string',{'$\sigma(W_{p}*S_{y})$','$\sigma(W_{\Delta}*T_{y})$','$\sigma(W_{p}*S_{y})+\sigma(W_{\Delta}*T_{y})$'});
xlabel('Frequency (rad/sec)','interpreter', 'tex');
ylabel('Gain (dB)','interpreter', 'tex');
grid
set(gca,'Color','w')

