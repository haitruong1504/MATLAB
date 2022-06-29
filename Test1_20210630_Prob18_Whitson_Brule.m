%pressure for 2PE calculation in Psia
%Temperature for 2PE calculation in oR
clc
format long
%% Input data
%Constant:
R = 10.7315; %psia*ft3/((lbm.mol)*R)

%Take example from Silva 2017
p = 500;    %psia
T = 740; %oR
%The order of zi is from C1, C4, C10
zi = [0.5; 0.42; 0.08];

%Input data of critical properties for components from Whitson and Brule
%The order is from H2O, CO2, C1, C2, C3, nC4, nC5, nC6, nC7 
Pc =[667.8; 550.7; 304.0]; %psia
Tc =[343.0; 765.3; 1111.8]; %oR
omega =[0.0115; 0.1928; 0.4902];

kij = [ 0,  0,  0;
        0,  0,  0;
        0,  0,  0];
    
tol = 1e-7;
maxiter = 50;
%Initial assumption of Ki
Ki = wilson(Pc, Tc, p, T, omega);

[Kiv, yiv, yil, Fv, Zv, Zl, ] = vapor_liquid_2Peq(Ki, p, T, Pc, Tc, zi, omega, kij, tol, maxiter);
