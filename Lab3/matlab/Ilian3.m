% LAB 3 - Ilian
clear all; clc; close all;


s=tf('s');
G=1e4*(s+2)/(s+3)/(s+100)^2;
Hinf(G);


%% 
Fsim=F; Gsim=G ;
% Edit parameters in macro.m
macro