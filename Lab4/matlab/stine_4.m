% Lab 4 - Stine
clear all; clc; close all;

% ----------------- 3.1.1 -----------------------
% --- G min phase ---
s = tf('s');
sysmp = minreal(minphase);
Gmp = minreal(sysmp.C*(s*eye(size(sysmp.A)) - sysmp.A)^(-1)*sysmp.B + sysmp.D);
Gmp0 = evalfr(Gmp, 0);
W1 = inv(Gmp0);
Gmp_tilde = Gmp * W1;
bode(Gmp_tilde); grid on;
% ----------------- 3.1.2 -----------------------
phi = pi / 3; % phase margin
wc = 0.1; % cross over frequency
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


% ----------------- 3.1.1 -----------------------
% --- G non min phase ---
sysNONmp = minreal(nonminphase);
G_NON = minreal(sysNONmp.C*(s*eye(size(sysNONmp.A)) - sysNONmp.A)^(-1)*sysNONmp.B + sysNONmp.D);
G_NON0 = evalfr(G_NON, 0);
W1_non = inv(G_NON0);
G_NON_tilde = G_NON * W1_non;
figure()
bode(G_NON_tilde); grid on;
% ----------------- 3.1.2 -----------------------
wc_NON = 0.02; % cross over frequency
% ---- for f1 -------
g12 = G_NON_tilde(1,2);
T12 = tan(phi - pi/2 - angle(evalfr(g12,i*wc_NON))) / wc_NON ;
l12 = g12 * (1 + 1/(s*T12));
K12 = 1 / abs(evalfr(l12, i*wc_NON));
f12 = minreal(K12 * (1 + 1/(s*T12)));
% ---- for f2 -------
g21 = G_NON_tilde(2,1);
T21 = tan(phi - pi/2 - angle(evalfr(g21,i*wc_NON))) / wc_NON ;
l21 = g21 * (1 + 1/(s*T21));
K21 = 1 / abs(evalfr(l21, i*wc_NON));
f21 = minreal(K21 * (1 + 1/(s*T21)));

F_NON_tilde = [0 f12 ; f21 0];
F_NON = W1_non * F_NON_tilde;
% ----------------- 3.1.3 -----------------------
S_NON = inv(eye(size(G_NON)) + G_NON * F_NON );
T_NON = S_NON * G_NON * F_NON;

figure()
bode(S_NON); hold on; grid on;
bode(T_NON);