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
wc_NON = 0.02;     % rad/s crossover frequency  (NON-minumum phase case!)
%---------------------------
%       FUNCTIONS:
T =@(G, w) tan(phi - pi/2 - angle(evalfr(G, 1i*w))) / w;
L =@(G, w) G*(1 + 1/(s*T(G,w)));
K =@(G, w) 1/abs(evalfr(L(G,w), 1i*w));
f =@(G, w) K(G, w)*(1 + 1/(s*T(G, w)));
%---------------------------
G_mp = G;
f1_mp = minreal(f(G_mp(1,1), wc_mp));
f2_mp = minreal(f(G_mp(2,2), wc_mp));

f1_NON = minreal(f(G_NON(1,2), wc_NON));
f2_NON = minreal(f(G_NON(2,1), wc_NON));

disp('------MINIMUM PHASE-------')
F_mp = [f1_mp, 0;0, f2_mp]
disp('-----NON-MINIMUM PHASE----')
F_NON = [0, f1_NON;f2_NON, 0]

L_mp = minreal(G_mp*F_mp);
L_NON = minreal(G_NON*F_NON);


disp('plot figure 3211: bode(L_mp) minimum phase')
figure(3211)
bode(L_mp); grid on; title('3.2.1(1): Bode(L_{mp}) minimum phase');


disp('plot figure 3212: bode(L1)')
figure(3212)
bode(L_NON); grid on; title('3.2.1(2): Bode(L_{NON}) non-minimum phase');


%% 3.2.2
disp(' ')
disp('<strong>---------------------------------------------------</strong>')
disp('<strong>Task 3.2.2</strong>')

S_NON = (eye(2) + L_NON)^(-1);
T_NON = S_NON*L_NON;

S_mp = (eye(2) + L_mp)^(-1);
T_mp = S_mp*L_mp;

RGA_mp_0 =  evalfr(G_mp,0).*(pinv(evalfr(G_mp,0)).');
RGA_mp_wc =  evalfr(G_mp,1i*wc_mp).*(pinv(evalfr(G_mp,1i*wc_mp)).');
disp('RGA G_{mp}(0): ')
disp(RGA_mp_0)
disp('RGA G_{mp}(iwc): ')
disp(RGA_mp_wc)

RGA_NON_0 = evalfr(G_NON,0).*(pinv(evalfr(G_NON,0)).');
RGA_NON_wc = evalfr(G_NON,1i*wc_NON).*(pinv(evalfr(G_NON,1i*wc_NON)).');
disp('RGA G_{NON}(0): ')
disp(RGA_NON_0)
disp('RGA G_{NON}(iwc): ')
disp(RGA_NON_wc)

% disp('MP')
% disp('u1 - y1')
% stepinfo(L_mp(1,1)/(1+L_mp(1,1))) % u1 - y1
% disp('u1 - y2')
% stepinfo(L_mp(1,2)/(1+L_mp(1,2))) % u1 - y2
% disp('u2 - y1')
% stepinfo(L_mp(2,1)/(1+L_mp(2,1))) % u2 - y1
% disp('u2 - y2')
% stepinfo(L_mp(2,2)/(1+L_mp(2,2))) % u2 - y2
% 
% disp('NON')
% stepinfo(L_NON(1,1)/(1+L_NON(1,1))) % u1 - y1
% stepinfo(L_NON(1,2)/(1+L_NON(1,2)))
% stepinfo(L_NON(2,1)/(1+L_NON(2,1)))
% stepinfo(L_NON(2,2)/(1+L_NON(2,2)))

G = G_NON;
F = F_NON;
% disp('Calls closedloop')
% closedloop
disp('CALLS close all')
close all;
