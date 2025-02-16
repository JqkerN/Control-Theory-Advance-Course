% Lab 4 - Ilian
clear all; clc; close all;

%% 3.1.1
disp('<strong> Exercise 3.1.1: </strong>')
s = tf('s');

%%%%%%%% MINPHASE %%%%%%%%%
sysmp = minreal(minphase);
G_MP = minreal(sysmp.C*(s*eye(size(sysmp.A)) - sysmp.A)^(-1)*sysmp.B + sysmp.D);
G_tilde_MP = minreal(G_MP*evalfr(G_MP, 0)^(-1));

disp('PLOTTING: 3.1.1 static decoupling, MINPHASE')
figure(311)
bode(G_tilde_MP); grid on;
title('3.1.1 static decoupling, MINPHASE')

%%%%%%% NON-MINPHASE %%%%%%%
sys_non_mp = minreal(nonminphase);
G_NON = minreal(sys_non_mp.C*(s*eye(size(sys_non_mp.A)) - sys_non_mp.A)^(-1)*sys_non_mp.B + sys_non_mp.D);
G_tilde_NON = minreal(G_NON*evalfr(G_NON, 0)^(-1));

disp('PLOTTING: 3.1.2 static decoupling, NON-MINPHASE')
figure(312)
bode(G_tilde_NON); grid on;
title('3.1.1 static decoupling, NON-MINPHASE')

%% 3.1.2
disp('<strong> Exercise 3.1.2: </strong>')

%---------------------------
%       VARIABLES:
phi = pi/3; % Intended phase margin
wc_MP = 0.1;   % rad/s crossover frequency  (minimum phase case!)
wc_NON = 0.02;     % rad/s crossover frequency  (NON-minumum phase case!)
%---------------------------
%       FUNCTIONS:
T =@(G, w) tan(phi - pi/2 - angle(evalfr(G, 1i*w))) / w;
L =@(G, w) G*(1 + 1/(s*T(G,w)));
K =@(G, w) 1/abs(evalfr(L(G,w), 1i*w));
f =@(G, w) K(G, w)*(1 + 1/(s*T(G, w)));
%---------------------------

VERBOSE = false;
if VERBOSE == true % TO IDENTIFY HOW TO MATCH F. 
    RGA_MP_0 =  evalfr(G_tilde_MP,0).*(pinv(evalfr(G_tilde_MP,0)).');
    RGA_MP_wc =  evalfr(G_tilde_MP,1i*wc_MP).*(pinv(evalfr(G_tilde_MP,1i*wc_MP)).');
    disp('RGA G_{mp}(0): ')
    disp(RGA_MP_0)
    disp('RGA G_{mp}(iwc): ')
    disp(RGA_MP_wc)

    RGA_NON_0 = evalfr(G_tilde_NON,0).*(pinv(evalfr(G_tilde_NON,0)).');
    RGA_NON_wc = evalfr(G_tilde_NON,1i*wc_NON).*(pinv(evalfr(G_tilde_NON,1i*wc_NON)).');
    disp('RGA G_{NON}(0): ')
    disp(RGA_NON_0)
    disp('RGA G_{NON}(iwc): ')
    disp(RGA_NON_wc)
end

f1_MP = minreal(f(G_tilde_MP(1,1), wc_MP));
f2_MP = minreal(f(G_tilde_MP(2,2), wc_MP));
f1_NON = minreal(f(G_tilde_NON(1,2), wc_NON));
f2_NON = minreal(f(G_tilde_NON(2,1), wc_NON));

F_tilde_MP = [f1_MP, 0;0, f2_MP];
F_MP = evalfr(G_MP, 0)^(-1)*F_tilde_MP;

F_tilde_NON = [0, f1_NON;f2_NON, 0];
F_NON = evalfr(G_NON, 0)^(-1)*F_tilde_NON;

if VERBOSE == true
    disp('------MINIMUM PHASE-------')
    disp(F_MP)
    disp('-----NON-MINIMUM PHASE----')
    disp(F_NON)
end

%% 3.1.3
disp('<strong> Exercise 3.1.3: </strong>')

S_MP = (eye(size(G_MP)) + G_MP*F_MP)^(-1); 
T_MP = (eye(size(G_MP)) + G_MP*F_MP)^(-1)*G_MP*F_MP; 
disp('PLOTTING: 3.1.3(1) Sensitivity and complementary, MINPHASE')
figure(3131)
bode(S_MP); hold on;
title('3.1.3(1) Sensitivity and complementary, MINPHASE'); grid on;
bode(T_MP); legend('S', 'T');

S_NON = (eye(size(G_NON)) + G_NON*F_NON)^(-1); 
T_NON = (eye(size(G_NON)) + G_NON*F_NON)^(-1)*G_NON*F_NON; 
disp('PLOTTING: 3.1.3(2) Sensitivity and complementary, NON-MINPHASE')
figure(3132)
bode(S_NON); hold on;
title('3.1.3(2) Sensitivity and complementary, NON-MINPHASE'); grid on;
bode(T_NON); legend('S', 'T');











