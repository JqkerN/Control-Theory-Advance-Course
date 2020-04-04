close all; clear all; clc;


%% 3.1.1
disp('<strong>Task 3.1.1</strong>')
s = tf('s');
sysmp = minreal(minphase);
G = minreal(sysmp.C*(s*eye(size(sysmp.A)) - sysmp.A)^(-1)*sysmp.B + sysmp.D);

size(G)
disp('Poles G(1,1):')
disp(pole(G(1,1)))
disp('Poles G(1,2):')
disp(pole(G(1,2)))
disp('Poles G(2,1):')
disp(pole(G(2,1)))
disp('Poles G(2,2):')
disp(pole(G(2,2)))
disp('Zeros G(1,1):')
disp(tzero(G(1,1)))
disp('Zeros G(1,2):')
disp(tzero(G(1,2)))
disp('Zeros G(2,1):')
disp(tzero(G(2,1)))
disp('Zeros G(2,2):')
disp(tzero(G(2,2)))

%% 3.1.2
disp(' ')
disp('<strong>---------------------------------------------------</strong>')
disp('<strong>Task 3.1.2</strong>')

disp('Poles G:')
disp(pole(G))
disp('Zeros G:')
disp(tzero(G))

%% 3.1.3
disp(' ')
disp('<strong>---------------------------------------------------</strong>')
disp('<strong>Task 3.1.3</strong>')

disp('plot figure 313: singular values')
figure(313)
sigma(G); title('3.1.3 singular values'); grid on;

%% 3.1.4
disp(' ')
disp('<strong>---------------------------------------------------</strong>')
disp('<strong>Task 3.1.4</strong>')

RGA = G.*inv(G)';
disp('RGA: ')
disp(evalfr(RGA,0))

%% 3.1.5
disp(' ')
disp('<strong>---------------------------------------------------</strong>')
disp('<strong>Task 3.1.5</strong>')

disp('plot figure 315: G(*,*) step response')
figure(315)
step(G); title('3.1.5 Step response'); grid on;
disp('Ans. One could notice that they are coupled. Which matches with the matrix from 3.1.5')

%% 3.1.6
disp(' ')
disp('<strong>---------------------------------------------------</strong>')
disp('<strong>Task 3.1.6</strong>')
s = tf('s');
sysNONmp = minreal(nonminphase);
G_NON = minreal(sysNONmp.C*(s*eye(size(sysNONmp.A)) - sysNONmp.A)^(-1)*sysNONmp.B + sysNONmp.D);

size(G_NON)
disp('Poles G_NON(1,1):')
disp(pole(G_NON(1,1)))
disp('Poles G_NON(1,2):')
disp(pole(G_NON(1,2)))
disp('Poles G_NON(2,1):')
disp(pole(G_NON(2,1)))
disp('Poles G_NON(2,2):')
disp(pole(G_NON(2,2)))
disp('Zeros G_NON(1,1):')
disp(tzero(G_NON(1,1)))
disp('Zeros G_NON(1,2):')
disp(tzero(G_NON(1,2)))
disp('Zeros G_NON(2,1):')
disp(tzero(G_NON(2,1)))
disp('Zeros G_NON(2,2):')
disp(tzero(G_NON(2,2)))


RGA_NON = G_NON.*inv(G_NON)';
disp('RGA: ')
disp(evalfr(RGA_NON,0))

disp('plot figure 3161: singular values')
figure(3161)
sigma(G_NON); title('3.1.6(1) singular values'); grid on;

disp('plot figure 3162: G(*,*) step response')
figure(3162)
step(G_NON); title('3.1.6(2) Step response'); grid on;

%% 3.2.1
clc; close all;
disp(' ')
disp('<strong>---------------------------------------------------</strong>')
disp('<strong>Task 3.2.1</strong>')
% tf, zoom, append, inv & feedback.
% G, G_NON
phi = pi/3;
wc_min = 0.1;
wc_nm = 0.02;

% ----------------- MIN PHASE -----------------------
% ---- for f1 -------
g11 = G(1,1);
T1 = tan(phi - pi + pi/2 - angle(evalfr(g11,i*wc_min))) / wc_min ;
l11 = g11 * (1 + 1/(s*T1));
K1 = 1 / abs(evalfr(l11, i*wc_min));
f1 = K1 * (1 + 1/(s*T1));
% ---- for f2 -------
g22 = G(2,2);
T2 = tan(phi - pi + pi/2 - angle(evalfr(g22,i*wc_min))) / wc_min ;
l22 = g22 * (1 + 1/(s*T2));
K2 = 1 / abs(evalfr(l22, i*wc_min));
f2 = K2 * (1 + 1/(s*T2));

F1_min = [f1 0 ; 0 f2];
F2_min = [0 f1 ; f2 0];

L1 = G * F1_min;
L2 = G * F2_min;
figure(1)
subplot(1,2,1)
bode(L1); grid on;
subplot(1,2,2)
bode(L2); grid on;
suptitle('Min phase')

% ----------------- NON MIN PHASE -----------------------
% ---- for f1 -------
g11 = G_NON(1,1);
T1 = tan(phi - pi + pi/2 - angle(evalfr(g11,i*wc_nm))) / wc_nm ;
l11 = g11 * (1 + 1/(s*T1));
K1 = 1 / abs(evalfr(l11, i*wc_nm));
f1 = K1 * (1 + 1/(s*T1));
% ---- for f2 -------
g22 = G_NON(2,2);
T2 = tan(phi - pi + pi/2 - angle(evalfr(g22,i*wc_nm))) / wc_nm ;
l22 = g22 * (1 + 1/(s*T2));
K2 = 1 / abs(evalfr(l22, i*wc_nm));
f2 = K2 * (1 + 1/(s*T2));

F1_nm = [f1 0 ; 0 f2];
F2_nm = [0 f1 ; f2 0];

L1 = G_NON * F1_nm;
L2 = G_NON * F2_nm;
figure(2)
subplot(1,2,1)
bode(L1); grid on;
subplot(1,2,2)
bode(L2); grid on;
suptitle('Non min phase')

%% 3.2.2
clc; close all;
disp(' ')
disp('<strong>---------------------------------------------------</strong>')
disp('<strong>Task 3.2.2</strong>')

% ----------------- ALL F -----------------------
F_all = {F1_min, F2_min, F1_nm, F2_nm};
F = F_all{3}
% ----------------- ALL G -----------------------
G_all = {G , G_NON};
G = G_all{2};
% ----- Sensativity & comp. Sensativity ---------
S = inv((eye(2) + G*F));
T1 = S * G*F;

% 3.2.3
close all;
disp(' ')
disp('<strong>---------------------------------------------------</strong>')
disp('<strong>Task 3.2.3</strong>')

% Skala om for G_non min??
closedloop