close all; clear all; clc;


%% 3.1.1
disp('<strong>Task 3.1.1</strong>')
s = tf('s');
sysmp = minreal(minphase);
Gmp = minreal(sysmp.C*(s*eye(size(sysmp.A)) - sysmp.A)^(-1)*sysmp.B + sysmp.D);

size(Gmp)
disp('Poles G(1,1):')
disp(pole(Gmp(1,1)))
disp('Poles G(1,2):')
disp(pole(Gmp(1,2)))
disp('Poles G(2,1):')
disp(pole(Gmp(2,1)))
disp('Poles G(2,2):')
disp(pole(Gmp(2,2)))
disp('Zeros G(1,1):')
disp(tzero(Gmp(1,1)))
disp('Zeros G(1,2):')
disp(tzero(Gmp(1,2)))
disp('Zeros G(2,1):')
disp(tzero(Gmp(2,1)))
disp('Zeros G(2,2):')
disp(tzero(Gmp(2,2)))

%% 3.1.2
disp(' ')
disp('<strong>---------------------------------------------------</strong>')
disp('<strong>Task 3.1.2</strong>')

disp('Poles G:')
disp(pole(Gmp))
disp('Zeros G:')
disp(tzero(Gmp))

%% 3.1.3
disp(' ')
disp('<strong>---------------------------------------------------</strong>')
disp('<strong>Task 3.1.3</strong>')

disp('plot figure 313: singular values')
figure(313)
sigma(Gmp); title('3.1.3 singular values'); grid on;

%% 3.1.4
disp(' ')
disp('<strong>---------------------------------------------------</strong>')
disp('<strong>Task 3.1.4</strong>')

RGA = Gmp.*inv(Gmp)';
disp('RGA: ')
disp(evalfr(RGA,0))

%% 3.1.5
disp(' ')
disp('<strong>---------------------------------------------------</strong>')
disp('<strong>Task 3.1.5</strong>')

disp('plot figure 315: G(*,*) step response')
figure(315)
step(Gmp); title('3.1.5 Step response'); grid on;
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
close all;
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
g11 = Gmp(1,1);
T11 = tan(phi - pi/2 - angle(evalfr(g11,i*wc_min))) / wc_min ;
l11 = g11 * (1 + 1/(s*T11));
K11 = 1 / abs(evalfr(l11, i*wc_min));
f11 = minreal(K11 * (1 + 1/(s*T11)));
% f11 = K11 * (1 + 1/(s*T11));
% ---- for f2 -------
g22 = Gmp(2,2);
T22 = tan(phi - pi/2 - angle(evalfr(g22,i*wc_min))) / wc_min ;
l22 = g22 * (1 + 1/(s*T22));
K22 = 1 / abs(evalfr(l22, i*wc_min));
f22 = minreal(K22 * (1 + 1/(s*T22)));
% f22 = K22 * (1 + 1/(s*T22));

F1_min = [f11 0 ; 0 f22];

L1 = Gmp * F1_min;
figure(1)
bode(L1); grid on; hold on;


% ----------------- NON MIN PHASE -----------------------
% ---- for f1 -------
g12 = G_NON(1,2);
T12 = tan(phi - pi/2 - angle(evalfr(g12,i*wc_nm))) / wc_nm ;
l12 = g12 * (1 + 1/(s*T12));
K12 = 1 / abs(evalfr(l12, i*wc_nm));
f12 = minreal(K12 * (1 + 1/(s*T12)));
% ---- for f2 -------
g21 = G_NON(2,1);
T21 = tan(phi - pi/2 - angle(evalfr(g21,i*wc_nm))) / wc_nm ;
l21 = g21 * (1 + 1/(s*T21));
K21 = 1 / abs(evalfr(l21, i*wc_nm));
f21 = minreal(K21 * (1 + 1/(s*T21)));

F2_nm = [0 f12 ; f21 0];

L2 = G_NON * F2_nm;

bode(L2); 
legend('Min phase' , 'Non min phase')
hold off;


%% 3.2.2
close all;
disp(' ')
disp('<strong>---------------------------------------------------</strong>')
disp('<strong>Task 3.2.2</strong>')

% ----------------- ALL F -----------------------
F_all = {F1_min, F2_nm};
F = F_all{1};
% ----------------- ALL G -----------------------
G_all = {Gmp , G_NON};
G = G_all{1};
% ----- Sensativity & comp. Sensativity ---------
S = inv((eye(2) + Gmp*F));
T = S * Gmp*F;

% 3.2.3
close all;
disp(' ')
disp('<strong>---------------------------------------------------</strong>')
disp('<strong>Task 3.2.3</strong>')

% closedloop