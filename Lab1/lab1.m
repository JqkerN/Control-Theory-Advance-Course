% LAB 1

% ------ CONTROL FUNCTIONS --------
% s = tf('s');
% % The product of two transfer functions is obtained by
% G12 = G1 * G2
% % For a system with 2 inputs and 2 outputs, the closed-loop transfer matrix is obtained with
% S=feedback(eye(2),G*F); T=feedback(G*F,eye(2))
% % For a SISO system this can be written
% S=1/(1+G*F); T=G*F/(1+G*F)
% % The bode diagram for G is plotted by typing
% bode(G) or bode(G,{wmin,wmax})
% % Amplitude and phase at a given frequency are obtained by
% [m,p]=bode(G,w)
% % Phase margin, amplitude margin and corresponding frequencies are obtained by
% [Gm,Pm,wp,wc]=margin(G*F)
% % To simulate a step response in the control signal, use the function
% step(G) or step(G,tfinal)
% % In the same way, to simulate a step response in the reference signal, we type
% step(T)
% -------------------------------

%% --- 4.1.1 ---
clear all; clc; close all;

s = tf('s');
G = 3*(-s + 1) / ((5*s+1)*(10*s+1));
% wc = 0.2; % [rad/s]
% phi = 30 * pi/180; % [rad]
% ---------- From Bode ----------
dB = -9.39;
mag = db2mag(dB);
phase = 30;% + 6 + 250-180 % gives TWO beta = 0.12
wc = 0.2;
wcd = 2*wc;
% --------------------------------

% ---------- Parameters ----------
beta = 0.1; %beta 4.1.1 = 0.33, beta 4.1.3 = 0.1;
gamma = 0;
td = 1 / (wcd*sqrt(beta));
ti = 10 / wcd;
% --------------------------------

% ---------- Controller ----------
K = sqrt(beta) / mag;
e1 = 1 / (3*K); % 3 when s->0 for G(s)
Lead = (td*s + 1) / (beta*td*s + 1);
Lag = (ti*s + 1) / (ti*s + gamma);
F = K * Lead * Lag;
% -------------------------------

% ---------- Plot ---------------
figure(411)
bode(G); hold on; grid on
bode(G*F); legend('G(s)' , 'G(s)*F(s)')
% --------------------------------


% --- 4.1.2 ---
% clear all; clc; close all;
figure(412)
Go = F*G;
Gc = Go / (1 + Go); 
step(Gc); grid on
stepinfo(Gc)

ymax = 1.23; yf= 1;
Mt = (ymax-yf)/yf;



% --- 4.1.3 ---
% clear all; clc; close all;


