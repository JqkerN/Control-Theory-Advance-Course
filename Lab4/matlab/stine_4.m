% Lab 4 - Stine
clear all; clc; close all;

% ----------------- 3.1.1 -----------------------
% --- G min phase ---
s = tf('s');
sysmp = minreal(minphase);
Gmp = minreal(sysmp.C*(s*eye(size(sysmp.A)) - sysmp.A)^(-1)*sysmp.B + sysmp.D);
Gmp0 = evalfr(Gmp, 0);
% W1 = inv(Gmp0);
phi = pi / 3; % phase margin
wc = 0.1; % cross over frequency

% --- 3.2.---
W1_12 = - Gmp(1,2)  / Gmp(1,1); %-0.0004446*(s + 0.05645)/0.0348/(s^2 + 0.07317*s + 0.001105);
W1_21 = - Gmp(2,1) / Gmp(2,2); %-0.0004649*(s + 0.05187)/0.03013/(s^2 + 0.08217*s + 0.001452);
W1 = [1, W1_12; W1_21 1];
W1 = minreal(W1);
% W1 = W1 * ((10*wc) / (s + 10*wc));
% --- 3.2 ---

Gmp_tilde = Gmp * W1;
bode(Gmp_tilde); grid on;
title('$\tilde{G}_{mp}$','Interpreter','Latex')
% ----------------- 3.1.2 -----------------------

% ---- for f1 -------
g11 = Gmp_tilde(1,1);
T11 = tan(phi - pi/2 - angle(evalfr(g11,i*wc))) / wc ;
l11 = g11 * (1 + 1/(s*T11));
K11 = 1 / abs(evalfr(l11, i*wc));
f11 = minreal(K11 * (1 + 1/(s*T11)));
% ---- for f2 -------
g22 = Gmp_tilde(2,2);
T22 = tan(phi - pi/2 - angle(evalfr(g22,i*wc))) / wc ;
l22 = g22 * (1 + 1/(s*T22));
K22 = 1 / abs(evalfr(l22, i*wc));
f22 = minreal(K22 * (1 + 1/(s*T22)));

Fmp_tilde = [f11 0 ; 0 f22];
Fmp = W1 * Fmp_tilde;
% ----------------- 3.1.3 -----------------------
Smp = inv(eye(size(Gmp)) + Gmp * Fmp );
Tmp = Smp * Gmp * Fmp;
figure()
bode(Smp); hold on; grid on;
bode(Tmp);
title('$S_{mp}$ and $T_{mp}$ ','Interpreter','Latex')





% ----------------- 3.1.1 -----------------------
% --- G non min phase ---
sysNONmp = minreal(nonminphase);
G_NON = minreal(sysNONmp.C*(s*eye(size(sysNONmp.A)) - sysNONmp.A)^(-1)*sysNONmp.B + sysNONmp.D);
G_NON0 = evalfr(G_NON, 0);
% W1_non = inv(G_NON0);
wc_NON = 0.02; % cross over frequency

% --- 3.2.---
W1_non_12 = - G_NON(1,2) / G_NON(1,1);%-0.003163*(s + 0.05106)/0.02088/(s^2 + 0.1378*s + 0.004265);
W1_non_21 = - G_NON(2,1) / G_NON(2,2); %-0.002586*(s + 0.04692)/0.01808/(s^2 + 0.1369*s + 0.004382);
% W1_non = [1 W1_non_12; W1_non_21 1];
W1_non_11 = - G_NON(2,2) / G_NON(2,1);
W1_non_22 = - G_NON(1,1) / G_NON(1,2);
W1_non = [W1_non_11 1;1 W1_non_22];
W1_non = W1_non * ((10*wc_NON) / (s + 10*wc_NON));
W1_non = minreal(W1_non);
% --- 3.2 ---

G_NON_tilde = G_NON * W1_non;
figure()
bode(G_NON_tilde); grid on;
title('$\tilde{G}_{NON}$','Interpreter','Latex')
% ----------------- 3.1.2 -----------------------
% ---- for f1 -------
g11 = G_NON_tilde(1,1);
T11 = tan(phi - pi/2 - angle(evalfr(g11,i*wc_NON))) / wc_NON ;
l11 = g11 * (1 + 1/(s*T11));
K11 = 1 / abs(evalfr(l11, i*wc_NON));
f11 = minreal(K11 * (1 + 1/(s*T11)));
% ---- for f2 -------
g22 = G_NON_tilde(2,2);
T22 = tan(phi - pi/2 - angle(evalfr(g22,i*wc_NON))) / wc_NON ;
l22 = g22 * (1 + 1/(s*T22));
K22 = 1 / abs(evalfr(l22, i*wc_NON));
f22 = minreal(K22 * (1 + 1/(s*T22)));

F_NON_tilde = [f11 0 ; 0 f22];
F_NON = W1_non * F_NON_tilde;
% ----------------- 3.1.3 -----------------------
S_NON = inv(eye(size(G_NON)) + G_NON * F_NON );
T_NON = S_NON * G_NON * F_NON;

figure()
bode(S_NON); hold on; grid on;
bode(T_NON);
title('$S_{NON}$  and $T_{NON}$ ','Interpreter','Latex')


%% ----------------- 3.1.4 -----------------------
G_all = {Gmp G_NON};
F_all = {Fmp F_NON};

G = G_all{1};
F = F_all{1};
closedloop


