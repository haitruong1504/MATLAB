%pressure for 2PE calculation in Psia
%Temperature for 2PE calculation in oR
clear;
clc;
format long
%% Input data
%Constant:
R = 10.7315; %psia*ft3/((lbm.mol)*R)

%Take example from Silva 2017
p = 1000;    %psia
T = 600; %oR
%The order of zi is from N2, CO2 C1, C2, C3, iC4, nC4, iC5, nC5, C6
zi = [0.001; 0.7895; 0.0344; 0.0085; 0.0063; 0.0008; 0.0069; 0.0033; 0.0045; 0.1448];

%Input data of critical properties for components from Whitson and Brule
%The order is from N2, CO2 C1, C2, C3, iC4, nC4, iC5, nC5, C6
Pc =[493.0; 1070.6; 667.8; 707.8; 616.3; 529.1; 550.7; 490.4; 488.6; 436.9]; %psia
Tc =[227.3; 547.6; 343.0; 549.8; 665.7; 734.7; 765.3; 828.8; 845.4; 913.4]; %oR
omega =[0.04; 0.225; 0.008; 0.098; 0.152; 0.176; 0.193; 0.227; 0.251; 0.296];

kij = [ 0,         -0.0315,     0.0278,     0.0407,     0.0763,     0.0944,     0.07,   0.0867,     0.0878,     0.08;
        -0.0315,    0,          0.12,       0.12,       0.12,       0.12,       0.12,   0.12,       0.12,       0.12;
        0.0278,     0.12,       0,          0,          0,          0,          0,      0,          0,          0;
        0.0407,     0.12,       0,          0,          0,          0,          0,      0,          0,          0;
        0.0763,     0.12,       0,          0,          0,          0,          0,      0,          0,          0;
        0.0944,     0.12,       0,          0,          0,          0,          0,      0,          0,          0;
        0.07,       0.12,       0,          0,          0,          0,          0,      0,          0,          0;
        0.0867,     0.12,       0,          0,          0,          0,          0,      0,          0,          0;
        0.0878,     0.12,       0,          0,          0,          0,          0,      0,          0,          0;
        0.08,       0.12,       0,          0,          0,          0,          0,      0,          0,          0;];
tol = 1e-6;
maxiter = 1e6;

[twophase, Kiv, Kil] = stability_test_on_single_phase(zi, p, T, Pc, Tc, omega, kij);

%[Kiv, yiv, yil, Fv, Zv, Zl] = vapor_liquid_2Peq(p, T, Pc, Tc, zi, omega, kij, tol, maxiter); %2 phase
