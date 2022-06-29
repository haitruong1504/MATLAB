%pressure for 2PE calculation in Psia
%Temperature for 2PE calculation in oR
% Tinh va so sanh voi PVTsim
format long;
clc;
clear;
%% Input data
%Constant:
R = 10.7315; %psia*ft3/((lbm.mol)*R)

%Take example from Silva 2017
p = 600;    %psia
T = 700; %oR
%The order of zi is from CO2 C1, C2, C3, nC4, nC5, C10
zi = [0.05; 0.1; 0.12; 0.12; 0.15; 0.17; 0.29];


%Input data of critical properties for components from Whitson and Brule
%The order is from CO2 C1, C2, C3, nC4, nC5, C10
Pc =[1069.87; 667.2; 708.35; 615.76; 551.1; 489.38; 322.73]; %psia
Tc =[547.56; 343.08; 549.72; 665.64; 765.36; 845.28; 1105.916]; %oR
omega =[0.225; 0.008; 0.098; 0.152; 0.193; 0.251; 0.4895];

kij = [ 0,      0.12,  0.12, 0.12, 0.12, 0.12,    0.1;
        0.12,   0,     0,    0,    0,    0,       0;
        0.12,   0,     0,    0,    0,    0,       0;
        0.12,   0,     0,    0,    0,    0,       0;
        0.12,   0,     0,    0,    0,    0,       0;
        0.12,   0,     0,    0,    0,    0,       0;
        0.1,    0,     0,    0,    0,    0,       0];
tol = 1e-8; % Fixed tolerance
maxiter = 50; % Fixed maxiter

[twophase, Ki_result, Yi] = stability_test_on_single_phase_20210828(zi, p, T, Pc, Tc, omega, kij);

% [phi_yi_1, ~] = fugacitycoef_multicomp(Yi(:,1), p, T, Pc, Tc, omega, kij);
% [phi_yi_2, ~] = fugacitycoef_multicomp(Yi(:,2), p, T, Pc, Tc, omega, kij);
% 
if twophase == true
    % For two phase.
    Ki_2 = Ki_result(:,1);
    [Kiv_2, yiv_2, yil_2, Fv_2, Zv_2, Zl_2] = vapor_liquid_2Peq(Ki_2, p, T, Pc, Tc, zi, omega, kij, tol, maxiter);
else
    % For one phase.
    [phi_overall, zfactor] = fugacitycoef_multicomp(zi, p, T, Pc, Tc, omega, kij);
end

% 
% First stability test on vapor phase composition
% [threephase_test_on_vapor, Yi_1_vapor, Yi_2_vapor] = stability_test_on_two_phase(yiv_2, yil_2, p, T, Pc, Tc, omega, kij);
% 
% % Second stability test on liquid phase composition
% [threephase_test_on_liquid, Yi_1_liquid, Yi_2_liquid] = stability_test_on_two_phase(yil_2, yiv_2, p, T, Pc, Tc, omega, kij);
% 
% if threephase_test_on_vapor == true || threephase_test_on_liquid == true
%     fprintf("Final conclusion is that the system has three phase\n");
% else
%     fprintf("Final conclusion is that the system has only two phase\n");
% end
% 
% Ki_3 = zeros(size(zi,1),2);
% 
% Case one of initial set
% Ki_3(:,1) = Kiv_2;
% Ki_3(:,2) = Yi_1_liquid./yil_2;
% 
% % Case two of initial set
% Ki_3(:,1) = Kiv_2;
% Ki_3(:,2) = Yi_2_liquid./yil_2;
% 
% Case three of initial set
% Ki_3(:,1) = Yi_1_liquid./yil_2;
% Ki_3(:,2) = Kiv_2;

% % Case four of initial set
% Ki_3(:,1) = Yi_2_liquid./yil_2;
% Ki_3(:,2) = Kiv_2;


% Case fifth of initial set
% Ki_3(:,1) = yiv_2./Yi_1_liquid;
% Ki_3(:,2) = yil_2./Yi_1_liquid;

% Case sixth of initial set
% Ki_3(:,1) = yiv_2./Yi_2_liquid;
% Ki_3(:,2) = yil_2./Yi_2_liquid;
% 
% % % Case seventh of initial set
% % Ki_3(:,1) = yil_2./Yi_1_liquid;
% % Ki_3(:,2) = yiv_2./Yi_1_liquid;
% 
% % Case seventh of initial set
% Ki_3(:,1) = yil_2./Yi_2_liquid;
% Ki_3(:,2) = yiv_2./Yi_2_liquid;
% 
% 
% [Ki_3, yiv_3, yiL2_3, yiL1_3, F_3] = vapor_liquid1_liquid2_3Peq(Ki_3, p, T, Pc, Tc, zi, omega, kij, tol, maxiter);







