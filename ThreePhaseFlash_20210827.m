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
p = 600; %psia
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

%% First stability test
fprintf("Stability test for one phases system. \n\n");
[twophases, Ki_first_stability_result, Yi_first_stability_result] = stability_test_on_single_phase_20210827(zi, p, T, Pc, Tc, omega, kij);

%% Two phase split calculations and second stability test
if twophases == true
    fprintf("\nSystem will split into two phases and need to perform two phase flash calculation.\n");
    % For two phase.
    Ki_2 = Ki_first_stability_result(:,1);
    [Kiv_2, yiv_2, yil_2, Fv_2, Zv_2, Zl_2] = vapor_liquid_2Peq(Ki_2, p, T, Pc, Tc, zi, omega, kij, tol, maxiter);
else
    fprintf("\nSystem is stable at single phase. \n");
    % For one phase.
    [phi_overall, zfactor] = fugacitycoef_multicomp(zi, p, T, Pc, Tc, omega, kij);
end

%% Second stability test
fprintf("\nStability test for two phases system. \n");
% % First stability test on vapor phase composition
fprintf("\nFirst stability test on vapor phase composition. \n\n");
[threephases_test_on_vapor, Ki_vapor_result, Yi_vapor_result] = stability_test_on_two_phase_20210827(yiv_2, p, T, Pc, Tc, omega, kij);

if threephases_test_on_vapor == false 
    fprintf("\nThe vapor phase of the two phases system is stable. \n");
else
    fprintf("\nThe vapor phase of the two phases system is not stable. \n");
end

% % Second stability test on liquid phase composition
fprintf("\nSecond stability test on vapor phase composition. \n\n");
[threephases_test_on_liquid, Ki_liquid_result, Yi_liquid_result] = stability_test_on_two_phase_20210827(yil_2, p, T, Pc, Tc, omega, kij);

if threephases_test_on_liquid == false
    fprintf("\nThe liquid phase of the two phases system is stable. \n");
else
    fprintf("\nThe liquid phase of the two phases system is not stable. \n");
end

%% Three phase spit calculation
Ki_3 = zeros(size(zi,1),2);

if threephases_test_on_vapor == true || threephases_test_on_liquid == true
    fprintf("\nSystem will split into three phases and need to perform three phase flash calculation\n");
    if threephases_test_on_vapor == true
        Ki_3(:,1) = Ki_vapor_result(:,1);
        Ki_3(:,2) = Kiv_2;
        [Ki_3, yiv_3, yiL2_3, yiL1_3, F_3] = vapor_liquid1_liquid2_3Peq(Ki_3, p, T, Pc, Tc, zi, omega, kij, tol, maxiter); 
    else
        Ki_3(:,1) = Ki_liquid_result(:,1);
        Ki_3(:,2) = Kiv_2;
        [Ki_3, yiv_3, yiL2_3, yiL1_3, F_3] = vapor_liquid1_liquid2_3Peq(Ki_3, p, T, Pc, Tc, zi, omega, kij, tol, maxiter);
    end
    
    if max(F_3) > 1
        fprintf("\nThe solution of three phase flash calculation is negative flash. Final conclusion is that the system has only two phases.\n");   
    else
        fprintf("\nFinal conclusion is that the system has three phases.\n");
    end
else
    fprintf("\n Final conclusion is that the system has only two phases.\n");
end
