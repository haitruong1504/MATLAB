%% Important notes
% Thay đổi initial guess cho stability test so với phiene bản 20210827
% Thay đổi yi và Yi trong file Check stability one phase

% Có cập nhật thêm ngày 29/8/2021

format long;
clc;
clear;
%% Input data
%Constant:
R = 10.7315; %psia*ft3/((lbm.mol)*R)

% Random example 
p1 = 500; %psia
T1 = 400; %oR

p2 = 1450.38; %psia
T2 = 635.67; %oR

p(1,1) = p1;
T(1,1) = T1;
for i = 2:20
    p(i,1) = (p2 - p1)/19 + p(i-1,1);
    T(i,1) = (T2 - T1)/19 + T(i-1,1);
end
% Thứ tự của CO2 luôn ở vị trí thứ 3

%The order of zi is from H2O, N2, CO2, C1, C2, C3, iC4, nC4, iC5, nC5, C6,
%C7, C8, C9, C10, C16
zi = [0.1; 0.00545; 0.02827; 0.45536; 0.08606; 0.0576; 0.01012; 0.02441; 0.00896; 0.01244; 0.01581; 0.02551; 0.02746; 0.01698; 0.01254; 0.11303];

%Input data of critical properties for components from PVTsim database
%The order is from H2O, N2, CO2, C1, C2, C3, iC4, nC4, iC5, nC5, C6,
%C7, C8, C9, C10, C16
Pc =[3203.73; 492.32; 1069.87; 667.2; 708.35; 615.76; 529.06; 551.10; 490.85; 489.38; 430.59; 427.16; 397.63; 363.54; 337.87; 248.91]; %psia
Tc =[1165.140; 227.16; 547.56; 343.08; 549.72; 665.64; 734.58; 765.36; 828.72; 845.28; 913.32; 968.679; 1004.477; 1047.74; 1084.527; 1289.501]; %oR
omega =[0.344; 0.04; 0.225; 0.008; 0.098; 0.152; 0.176; 0.193; 0.227; 0.251; 0.296; 0.3374; 0.3743; 0.4205; 0.4628; 0.7124];
critical_volume = [0.90; 1.44; 1.51; 1.59; 2.37; 3.25; 4.21; 4.08; 4.90; 4.87; 5.93; 7.62; 7.81; 8.48; 9.19; 14.76];
MW = [18.015; 28.014; 44.010; 16.043; 30.070; 44.097; 58.124; 58.124; 72.151; 72.151; 86.178; 96.000; 107.000; 121.000; 134.000; 222.000];

% Binary interaction parameters is read from excel file with the order of
% components is H2O, N2, CO2, C1, C2, C3, iC4, nC4, iC5, nC5, C6,
% C7, C8, C9, C10, C16
% This binary interaction parameters is from Winprop manual
kij = xlsread('Kij.xlsx');

tol = 1e-8; % Fixed tolerance
maxiter = 50; % Fixed maxiter

for i = 1:20
[Ki_final_result{i}, yi_final_result{i}, F_final_result{i}] = threephase_flash_calculation(zi, p(i,1), T(i,1), Pc, Tc, omega, MW, kij, tol, maxiter);
end
