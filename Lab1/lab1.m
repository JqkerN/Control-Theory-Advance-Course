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
bode(Gc); legend('G(s)' , 'Gc')
% --------------------------------


% --- 4.1.2 ---
% clear all; clc; close all;
figure(412)
step(Gc); grid on
stepinfo(Gc)

ymax = 1.23; yf= 1;
Mt = (ymax-yf)/yf;

%% --- 4.2.1 ---
clear all; clc; close all;

s = tf('s');
G = 20 / ( (s+1) * ( (s/20)^2 + s/20 + 1) );
Gd = 10 / (s+1);
wc = 10;
L = wc / s;
Fy = L*G^(-1); % Not proper

p1 = 10*wc;
p2 = 10*wc;
pole1 = 1/(s/p1 + 1);
pole2 = 1/(s/p2 + 1);
Fy_prop = Fy*pole1*pole2; % Proper

Go = G*Fy;
Gc = Go/(1 + Go);
Go_prop = G*Fy_prop;
Gc_prop = Go_prop/(1 + Go_prop);

% figure(4211)
% bode(Gc); hold on;
% bode(Gc_prop); legend('Gc not proper', 'Gc proper');
% figure(4212)
% step(Gc); hold on;
% step(Gc_prop); legend('Gc not proper', 'Gc proper');
% figure()
% bode(G); hold on; grid on;
% bode(Gd); legend('G' , 'Gd');
% figure()
% step(G); hold on; grid on;
% step(Gd); legend('G' , 'Gd');
% figure()
% pzmap(G); grid on;


% --- 4.2.2 ---
%--- Parameters -------
p1 = 50;
p2 = 50;
wI = 5;
Gd = [Gd];
% Gd2 = 1;
% Gd3 = G;
% --------------------

pole1 = 1/(s/p1 + 1);
pole2 = 1/(s/p2 + 1);

for ii = 1:1
    Fy1 = (s+wI)/s * G^(-1) * Gd(ii); % Unmodified
    Fy2 = Fy1*pole1*pole2; % Proper modified

%     Go = G*Fy1 + Gd(ii);
%     Gc1 = Go/(1 + Go);
%     Go = G*Fy2 + Gd(ii);
%     Gc2 = Go/(1 + Go);
    y1 = 1 / (1+G*Fy1)*Gd(ii);
    y2 = 1 / (1+G*Fy2)*Gd(ii);
%     subplot(2,2,ii)
% %     figure()
%     bode(Gc1); hold on;
%     bode(Gc2); legend('y unmodified', 'y modified');
%     subplot(1,2,ii)
% %     figure()
%     step(y1); hold on;
%     step(y2); legend('Gc', 'Gc modified');
%     stepinfo(y2)
    % figure(4223)
    % pzmap(Gc2); grid on;
end


% --- 4.2.3 ---
F = Fy2;
S = 1/ (1+L);
Go = F*G;
Gc = Go / (1 + Go); 
figure(1)
% subplot(1,2,1)
margin(Gc); hold on;
figure(2)
% subplot(1,2,2)
step(Gc); hold on;
% 1) decide wc and phase
wc = 14.3; phi = 30;
wcd = 20;
% phase = phi + 6 - 180
dB = -3.77;
mag = db2mag(dB);
%

% r = y+1;
% u = Fr*Fy*S - Fy2
% ---------- Parameters ----------
beta = 0.33; 
gamma = 0;
td = 1 / (wcd*sqrt(beta));
ti = 10 / wcd;
% --------------------------------

% ---------- Controller ----------
K = 1/(sqrt(beta)* mag);
e1 = 1 / (3*K); % 3 when s->0 for G(s)
Lead = (td*s + 1) / (beta*td*s + 1);
Lag = (ti*s + 1) / (ti*s + gamma);
tao = 0.01;
Fr = 1 / (1+ tao*s);
F = F * K * Lead;% * Lag;
Gc = G*Fr / (1 + G*Fy); 
figure(1)
% subplot(1,2,1)
margin(Gc); legend('Gc(s)', 'Gc(s) with Lead'); hold off;
figure(2)
% subplot(1,2,2)
step(Gc); legend('Gc(s)', 'Gc(s) with Lead'); hold off;
stepinfo(Gc)
% -------------------

