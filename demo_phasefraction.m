%% Phase mole fraction calculation
% The method used here is based on
% Okuno et al., 2010, A New Algorithm for Rachford-Rice For Multiphase
% Compositional Simulation, SPE Journal, 14 (2): 313-325

format long;
clear;
clc;
% % First example
% % Overall composition, 7 components.
% comp_overall = [
%     0.204322077
%     0.070970999
%     0.267194323
%     0.296291965
%     0.067046081
%     0.062489248
%     0.031685307 ];
% % Equilibrium constant, 3 phases.
% K = [
%     1.234669887  1.527133414
%     0.897277011  0.024564880
%     2.295257081  1.463482405
%     1.589548999  1.160905462
%     0.233493486  0.241662899
%     0.020381086  0.148152826
%     1.407156410  14.31280108 ];
% 
% % tolerance and the maximu iteration.
% tol = 1e-8;
% maxiter = 100;
% 
% [phasefrac, comp] = phasefraction(K, comp_overall, tol, maxiter);
% 
% disp("phasefrac is:"), disp(phasefrac);
% disp("comp is:"), disp(comp);

% 
% % Example 2:
% % Overall composition, 7 components.
% comp_overall = [
%     0.132266176697
%     0.205357472415
%     0.170087543100
%     0.186151796211
%     0.111333894738
%     0.034955417168
%     0.159847699672 ];
% % Equilibrium constant, 3 phases.
% K = [
%     26.3059904941  66.7435876079
%     1.91580344867  1.26478653025
%     1.42153325608  0.94711004430
%     3.21966622946  3.94954222664
%     0.22093634359  0.35954341233
%     0.01039336513  0.09327536295
%     19.4239894458  12.0162990083];
% 
% % tolerance and the maximu iteration.
% tol = 1e-8;
% maxiter = 100;
% 
% [phasefrac, comp] = phasefraction(K, comp_overall, tol, maxiter);
% 
% disp("phasefrac is:"), disp(phasefrac);
% disp("comp is:"), disp(comp);

% 
% % Example 3:
% % Overall composition, 7 components.
% comp_overall = [
%     0.896646630194
%     0.046757914522
%     0.000021572890
%     0.000026632729
%     0.016499094171
%     0.025646758089
%     0.014401397406 ];
% % Equilibrium constant, 3 phases.
% K = [
%     1.64571122126  1.61947897153
%     1.91627717926  2.65352105653
%     0.71408616431  0.68719907526
%     0.28582415424  0.18483049029
%     0.04917567928  0.01228448216
%     0.00326226927  0.00023212526
%     0.00000570946  0.00000003964];
% 
% % tolerance and the maximu iteration.
% tol = 1e-8;
% maxiter = 100;
% 
% [phasefrac, comp] = phasefraction(K, comp_overall, tol, maxiter);
% 
% disp("phasefrac is:"), disp(phasefrac);
% disp("comp is:"), disp(comp);


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


% Bao cao trong Báo cáo ngày 02 tháng 08 năm 2021
% 
% % Example 5: Table 2 from paper: Initialization of phase fractions in Rachford–Rice equations for robust and efficient three-phase split calculation
% % Overall composition, 3 components.
% comp_overall = [
%     0.66731
%     0.09575 
%     0.03540
%     0.04452 
%     0.08589
%     0.04470
%     0.02643 ];
% % Equilibrium constant, 3 phases.
% K = [
%     1.40089114681102        1.42336619958799
%     2.41359153035331        1.56360101270076
%     0.684675481993755       0.805778846552492
%     0.192706323169157       0.437918929556065
%     1.344316808771735e-2    0.136423337258229
%     2.913379631601974e-4    2.241151325196582e-2
%     9.614643893437818e-8    3.114699395928320e-4];
% 
% % tolerance and the maximu iteration.
% tol = 1e-8;
% maxiter = 100;
% 
% [phasefrac, comp] = phasefraction(K, comp_overall, tol, maxiter);
% 
% disp("phasefrac is:"), disp(phasefrac);
% disp("comp is:"), disp(comp);

% 
% % Example 6: Table 1 from paper: Initialization of phase fractions in Rachford–Rice equations for robust and efficient three-phase split calculation
% % Overall composition, 3 components.
% comp_overall = [
%     0.47 
%     0.126754033873246 
%     0.123759275876241 
%     0.190491864809508 
%     5.352678894647322e-2 
%     3.546803696453197e-2];
% % Equilibrium constant, 3 phases.
% K = [
%     0.886975564280731       1.85133355509695
%     183.729456216368        0.567851997436811
%     28.8439229979536        0.291644844783998
%     0.762796901964099       0.182989507250403 
%     6.805250689498878e-2    8.745408265736165e-2
%     0.345376016039736       0.623957189693138];
% 
% % tolerance and the maximu iteration.
% tol = 1e-8;
% maxiter = 100;
% 
% [phasefrac, comp] = phasefraction(K, comp_overall, tol, maxiter);
% 
% disp("phasefrac is:"), disp(phasefrac);
% disp("comp is:"), disp(comp);


% Example 7: Table 1 from paper: Initialization of phase fractions in Rachford–Rice equations for robust and efficient three-phase split calculation
% Overall composition, 3 components.
comp_overall = [
    0.466 
    0.127710667872289  
    0.124693307875307 
    0.191929538808070 
    5.393076494606923e-2 
    3.573572096426427e-2 ];
% Equilibrium constant, 3 phases.
K = [
    0.367489928755904       1.45983188593810
    91.9551101941298        0.627700554178016
    17.6437660816506        0.405472131110146
    0.523968443113866       0.291902855037650 
    5.444380423358842e-2    0.172272959622522
    0.192716832533260           0.704057279260822];

% tolerance and the maximu iteration.
tol = 1e-8;
maxiter = 100;

[phasefrac, comp] = phasefraction(K, comp_overall, tol, maxiter);

disp("phasefrac is:"), disp(phasefrac);
disp("comp is:"), disp(comp);