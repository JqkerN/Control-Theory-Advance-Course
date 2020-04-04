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
disp(' ')
disp('<strong>---------------------------------------------------</strong>')
disp('<strong>Task 3.2.1</strong>')

%---------------------------
%       VARIABLES:
phi = pi/3; % Intended phase margin
wc_mp = 0.1;   % rad/s crossover frequency  (minimum phase case!)
wc = 0.02;     % rad/s crossover frequency  (NON-minumum phase case!)
%---------------------------
%       FUNCTIONS:
T =@(G, w) tan(phi - pi/2 - angle(evalfr(G, 1i*w))) / w;
L =@(G, w) G*(1 + 1/(s*T(G,w)));
K =@(G, w) 1/abs(evalfr(L(G,w), 1i*w));
f =@(G, w) K(G, w)*(1 + 1/(s*T(G, w)));
%---------------------------
G_mp = G;
f1_mp = f(G_mp(1,1), wc_mp);
f2_mp = f(G_mp(2,2), wc_mp);

f1_NON = f(G_NON(1,1), wc);
f2_NON = f(G_NON(2,2), wc);

F_1_mp = [f1_mp, 0;0, f2_mp];
F_2_mp = [0, f1_mp;f2_mp, 0];

F_1_NON = [f1_NON, 0;0, f2_NON];
F_2_NON = [0, f1_NON;f2_NON, 0];

L1_mp = G_mp*F_1_mp;
L2_mp = G_mp*F_2_mp;
L1 = G_NON*F_1_NON;
L2 = G_NON*F_2_NON;

disp('plot figure 3211: bode(L1_mp) minphase')
figure(3211)
bode(L1_mp); grid on; title('3.2.1(1) Bode(L1_{mp}) minphase');
disp('plot figure 3212: bode(L2_mp) minphase')
figure(3212)
bode(L2_mp); grid on; title('3.2.1(2) Bode(L2_{mp}) minphase');

disp('plot figure 3213: bode(L1)')
figure(3213)
bode(L1); grid on; title('3.2.1(3) Bode(L1)');
disp('plot figure 3214: bode(L2)')
figure(3214)
bode(L2); grid on; title('3.2.1(4) Bode(L2)');


%% 3.2.2
disp(' ')
disp('<strong>---------------------------------------------------</strong>')
disp('<strong>Task 3.2.2</strong>')

S1 = (eye(2) + L1)^(-1);
S2 = (eye(2) + L2)^(-1);
T1 = S1*L1;
T2 = S2*L2;

S1_mp = (eye(2) + L1_mp)^(-1);
S2_mp = (eye(2) + L2_mp)^(-1);
T1_mp = S1_mp*L1_mp;
T2_mp = S2_mp*L2_mp;

G = G_NON;
F = F_1_NON;
