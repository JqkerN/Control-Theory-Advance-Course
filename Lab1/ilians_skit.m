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
Go = F*G;
Gc = Go / (1 + Go); 

% -------------------------------
% ---------- Plot ---------------
figure(411)
bode(G); hold on; grid on
bode(Gc); legend('G(s)' , 'Gc(s)')

% --------------------------------
% --- 4.1.2 ---
% clear all; clc; close all;
figure(412)
step(Gc); grid on
stepinfo(Gc)

ymax = 1.23; yf= 1;
Mt = (ymax-yf)/yf;

% --- 4.1.3 ---
% clear all; clc; close all;


%% ------------------------------
% ------------4.2----------------
% -------------------------------
s = tf('s');
G = 20/((s+1)*( (s/20)^2 + (s/20) +1));
Gd = 10/(s+1);

figure(42)
subplot(1,2,1)
bode(G); hold on;
bode(Gd); legend('G', 'Gd');

subplot(1,2,2)
step(G); hold on;
step(Gd); legend('G', 'Gd');

disp('Poles (G):')
disp(pole(G))
disp('Poles (Gd):')
disp(pole(Gd))

figure(1000)
pzmap(G)
grid on

%% ------------4.2.1--------------
% -------------------------------
clear all; clc; 
s = tf('s');
G = 20/((s+1)*( (s/20)^2 + (s/20) +1));
Gd = 10/(s+1);

wc = 10;
p1 = 10*wc;
p2 = 10*wc;

Fy = G^(-1)*10/s; % Not-proper
pole1 = 1/(s/p1 + 1);
pole2 = 1/(s/p2 + 1);
Fy_prop = Fy*pole1*pole2; % Proper

Go = G*Fy;
Gc = Go/(1 + Go);
Go_prop = G*Fy_prop;
Gc_prop = Go_prop/(1 + Go_prop);

figure(4211)
bode(Gc); hold on;
bode(Gc_prop); legend('Gc not proper', 'Gc proper');
figure(4212)
step(Gc); hold on;
step(Gc_prop); legend('Gc not proper', 'Gc proper');

disp('Gc:')
disp(stepinfo(Gc))

disp('Gc proper:')
disp(stepinfo(Gc_prop))



%% ------------4.2.2--------------
% -------------------------------
clear all; clc;

s = tf('s');
G = 20/((s+1)*( (s/20)^2 + (s/20) +1));
Gd = 10/(s+1);
w_I = 5;
wc = 10;

p1 = 50;
p2 = 500;
p3 = 1500;

pole1 = 1/(s/p1 + 1);
pole2 = 1/(s/p2 + 1);
pole3 = 1/(s/p3 + 1);


Fy_not_prop = (s + w_I)/s * G^(-1) * Gd;
Fy_prop_1 = Fy_not_prop*pole1*pole1 % Proper
Fy_prop_2 = Fy_not_prop*pole2*pole2; % Proper
Fy_prop_3 = Fy_not_prop*pole1*pole2*pole3; % Proper


Gc_not_prop = Gd/(1 + G*Fy_not_prop);
Gc_prop_1 = Gd/(1 + G*Fy_prop_1);
Gc_prop_2 = Gd/(1 + G*Fy_prop_2);
Gc_prop_3 = Gd/(1 + G*Fy_prop_3);



figure(4222)
step(minreal(Gc_not_prop)); hold on; grid on;
step(minreal(Gc_prop_1)); step(minreal(Gc_prop_2)); step(minreal(Gc_prop_3)); legend('Improper', 'Proper, p1=p2=50', 'Proper, p1=p2=500', 'Proper, p1=50 , p2=500, p3=1500');




%% ------------4.2.3--------------
% -------------------------------
clear all; clc;
p1 = 50;
s = tf('s');
G = 20/((s+1)*( (s/20)^2 + (s/20) +1));
Gd = 10/(s+1);
w_I = 5;
wc = 14.3;

pole1 = 1/(s/p1 + 1);
Fy_not_prop = (s + w_I)/s * G^(-1) * Gd;
Fy_prop_1 = Fy_not_prop*pole1*pole1; % Proper
Go = Fy_prop_1*G;
Gc = Go/(1+Go);

% Lead link
wcd = 20;
beta = 0.35;
tao_d = 1/(wcd*sqrt(beta));
K = 1/db2mag(-2);%/(sqrt(beta)*db2mag(-3.77));% sqrt(beta) / db2mag(-16);

Fy_lead = Fy_prop_1 * K* (tao_d*s + 1)/(beta*tao_d*s + 1);
Go_lead = Fy_lead*G;
Gc_lead = Go_lead/(1+Go_lead);

figure(4231)
subplot(1,2,1); bode(minreal(Gc)); grid on; hold on;
bode(minreal(Gc_lead)); 
subplot(1,2,2); step(minreal(Gc)); grid on; hold on;
step(minreal(Gc_lead));

disp('Gc:')
disp(stepinfo(minreal(Gc)))
disp('Gc_lead:')
disp(stepinfo(minreal(Gc_lead)))
%% Fr
tao = 0.09;
Fr = 1/(1+tao*s);
Gc_tot = Gd/(1+G*Fy_lead);

% Something
figure(4232)
step(minreal(Gc_tot));grid on;


% 
Gc_tot = G*Fy_lead*Fr/(1+G*Fy_lead);
figure(4231)
subplot(1,2,1); bode(minreal(Gc)); grid on; hold on;
bode(minreal(Gc_lead)); bode(minreal(Gc_tot));
subplot(1,2,2); step(minreal(Gc)); grid on; hold on;
step(minreal(Gc_lead));step(minreal(Gc_tot)); legend('Gc', 'Gc with lead', 'Gc with lead and Fr')
% 
% 
% bode(minreal(Gc_tot)); legend('Gc', 'Gc with lead', 'Gc with lead and Fr')
% subplot(1,2,2);
% step(minreal(Gc_tot)); legend('Gc', 'Gc with lead', 'Gc with lead and Fr')

disp('With F_lead and F_r:')
disp(stepinfo(minreal(Gc_tot)))

%% Size of control signal 
figure(4234)
step(Fy_lead*Fr/(1+Fy_lead*G)); hold on;
step(Fy_lead*Gd/(1+Fy_lead*G));
step(Fy_lead*Fr/(1+Fy_lead*G) - Fy_lead*Gd/(1+Fy_lead*G));
legend('FyFrSr','FyGdSd', 'FyFrSr - FyGdSd')


%% 4.2.4

figure(424)
bode(minreal(1/(1 + Fy_lead*G))); hold on;
bode(minreal(G*Fy_lead/(1+G*Fy_lead))); legend('Sensitivity', 'Complementary sensitivity');