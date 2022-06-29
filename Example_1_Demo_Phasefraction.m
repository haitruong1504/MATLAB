%Demo phase fraction with example from Efficient and Robust Three-Phase Split
%Computations (Abbas Firoozabadi 2010)

format long;
clear;
clc;

%First example     %FAIL FAIL FAIL FAIL FAIL
% Overall composition, 4 components.
comp_overall = [
    0.2;
    0.5;
    0.2;
    0.1;];
% Equilibrium constant, 3 phases.
K = [4.3449,  0.0003964;
    1.1366,  0.0001265;
    0.3042,  0.0002398;
    11.577,  8699.3;];

% tolerance and the maximum iteration.
tol = 1e-6;
maxiter = 100;

[phasefrac, comp] = phasefraction(K, comp_overall, tol, maxiter);

disp("phasefrac is:"), disp(phasefrac);
disp("comp is:"), disp(comp);

% % %Second example
% % % Overall composition, 5 components.
% % comp_overall = [
% %     0.66;
% %     0.03;
% %     0.01;
% %     0.05;
% %     0.25;];
% % % Equilibrium constant, 3 phases.
% % K = [9.7474,  6.9448;
% %     0.2834,  4.3275;
% %     0.0388,  5.8241;
% %     0.1776,  1.0070;
% %     0.0079,  0.1740];
% % 
% % % tolerance and the maximum iteration.
% % tol = 1e-6;
% % maxiter = 20;
% % 
% % [phasefrac, comp] = phasefraction(K, comp_overall, tol, maxiter);
% % 
% % disp("phasefrac is:"), disp(phasefrac);
% % disp("comp is:"), disp(comp);
% % 
% % 
% % %Third example
% % % Overall composition, 5 components.
% % comp_overall = [
% %     0.66;
% %     0.03;
% %     0.01;
% %     0.05;
% %     0.25;];
% % % Equilibrium constant, 3 phases.
% % K = [8.0930,  0;
% %     4.5132,   0;
% %     5.9477,   0;
% %     0.9646,   15.318;
% %     0.1329,  1.1972];
% % 
% % % tolerance and the maximum iteration.
% % tol = 1e-6;
% % maxiter = 100;
% % 
% % [phasefrac, comp] = phasefraction(K, comp_overall, tol, maxiter);
% % 
% % disp("phasefrac is:"), disp(phasefrac);
% % disp("comp is:"), disp(comp);
% 
% 
% % Example 4:
% % Overall composition, 3 components.
% comp_overall = [
%     0.08860
%     0.81514
%     0.09626 ];
% % Equilibrium constant, 3 phases.
% K = [
%     0.112359551  1.011235955
%     13.72549020  0.980392157
%     3.389830508  0.847457627];
% 
% % tolerance and the maximu iteration.
% tol = 1e-8;
% maxiter = 100;
% 
% [phasefrac, comp] = phasefraction(K, comp_overall, tol, maxiter);
% 
% disp("phasefrac is:"), disp(phasefrac);
% disp("comp is:"), disp(comp);
% Ki = [1;2;3];
% zi = [2;3;4];
% ncomp = size(zi,1);
% Ki_test (:,1) = Ki;
% Ki_test (:,2) = 1./Ki;
% Ki_test (:,3) = (Ki).^(1/3);
% Ki_test (:,4) = 1./((Ki).^(1/3));
% i = 5;
%     for j = 1:ncomp
%         for k = 1:ncomp
%             if k == j
%                 Ki_test(k,i) = (0.9)/zi(k,1);
%             else
%                 Ki_test(k,i) = (0.1/(ncomp-1))/zi(k,1);
%             end
%         end
%         i = i +1;
%     end
% zi = [1;2;3;4;5;6;7;8;9];
% ncomp = size(zi,1);
% 
% % Ki_test(3,5) = (0.9)/zi(3,1);
% 
% for i = 1: ncomp
%     if i == 3
%         Ki_test(i,1) = (0.9)/zi(i,1);
%     else
%         Ki_test(i,1) = (0.1/(ncomp-1))/zi(i,1);
%     end
% end


% 
% zi = [2;3;4];
% phi_zi = [1;2;3];
% 
% phi_Yi = [1, 2; 2, 3; 3, 4];
% 
% Yi_result = [3,4; 4,5; 5,6];
% 
% 
% TPD = [];
% % 
% % for i = 1: size(phi_Yi,2)
% %     TPD(i,1) = 1;
% %     for j = 1: size(zi,1)
% %         TPD(i,1) = TPD(i,1) + Yi_result(j,i)*(log(phi_Yi(j,i)) + log(Yi_result(j,i)) - log(phi_zi(j,1)) - log(zi(j,1)) - 1);
% %     end
% % end
% 
% 
% for i = 1: size(phi_Yi,2)
%     TPD(i,1) = 0;
%     for j = 1: size(zi,1)
%         TPD(i,1) = TPD(i,1) + Yi_result(j,i)*(log(phi_Yi(j,i)) + log(Yi_result(j,i)) - log(phi_zi(j,1)) - log(zi(j,1)));
%     end
% end
% 
% [~, index] = min(TPD);
