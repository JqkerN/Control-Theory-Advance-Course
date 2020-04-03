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
s = 0;
sysmp_zero = minreal(minphase);
G_RGA = sysmp_zero.C*(s*eye(size(sysmp_zero.A)) - sysmp_zero.A)^(-1)*sysmp_zero.B + sysmp_zero.D;
RGA_zero = G_RGA.*inv(G_RGA)'

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

s = 0;
sysNONmp_zero = minreal(nonminphase);
G_RGA_NON = sysNONmp_zero.C*(s*eye(size(sysNONmp_zero.A)) - sysNONmp_zero.A)^(-1)*sysNONmp_zero.B + sysNONmp_zero.D;
RGA_zero_NON = G_RGA_NON.*inv(G_RGA_NON)'

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