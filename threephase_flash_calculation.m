function [Ki_final_result, yi_final_result, F_final_result] = threephase_flash_calculation(zi, p, T, Pc, Tc, omega, MW, kij, tol, maxiter)

%% First stability test
fprintf("Stability test for one phases system. \n\n");
[twophases, Ki_first_stability_result, ~] = stability_test_on_single_phase_20210828(zi, p, T, Pc, Tc, omega, kij);

if twophases == true
    %% Two phase split calculations and second stability test
    fprintf("\nSystem will split into two phases and need to perform two phase flash calculation.\n");
    % For two phase.
    Ki_2 = Ki_first_stability_result;
    [twophase_flash, Kiv_2, yiv_2, yil_2, Fv_2] = vapor_liquid_2Peq(Ki_2, p, T, Pc, Tc, zi, omega, kij, tol, maxiter);
    Fl_2 = 1 - Fv_2;
    
    if twophase_flash == true
        %% Calculate phase molecular weight of each phase from two phase flash
        MW_vapor_2pflash = sum(yiv_2.*MW);
        MW_liquid_2pflash = sum(yil_2.*MW);
        
        %% Second stability test
        fprintf("\nStability test for two phases system. \n");
        if MW_vapor_2pflash > MW_liquid_2pflash
            fprintf("\nTwo phase stability test is conducted on vapor phase composition. \n\n");
            [threephases_test_on_vapor, Ki_vapor_result, ~] = stability_test_on_two_phase_20210828(yiv_2, Ki_first_stability_result, Kiv_2, p, T, Pc, Tc, omega, kij);
    
            if threephases_test_on_vapor == false 
                fprintf("\nThe two phases system is stable. \n");
                fprintf("\n Final conclusion is that the system has only two phases.\n");
        
                Ki_final_result = Kiv_2;
                yi_final_result(:,1) = yiv_2;
                yi_final_result(:,2) = yil_2;
                F_final_result(1,1) = Fv_2;
                F_final_result(2,1) = Fl_2;
        
            else
                fprintf("\nThe two phases system is not stable. \n");
                fprintf("\nSystem will split into three phases and need to perform three phase flash calculation\n");
        
                Ki_3(:,1) = 1./Kiv_2; % Two phase flash calculation
                Ki_3(:,2) = Ki_vapor_result; % Two phase stability test
                [threephase_flash, Ki_3, yiv_3, yiL2_3, yiL1_3, F_3] = vapor_liquid1_liquid2_3Peq(Ki_3, p, T, Pc, Tc, zi, omega, kij, tol, maxiter);
                Fl_3 = 1 - sum(F_3);
                
                if threephase_flash == false
                    fprintf("\n Final conclusion is that the system has only two phases.\n");
                    Ki_final_result = Kiv_2;
                    yi_final_result(:,1) = yiv_2;
                    yi_final_result(:,2) = yil_2;
                    F_final_result(1,1) = Fv_2;
                    F_final_result(2,1) = Fl_2;

                else
                    fprintf("\n Final conclusion is that the system has three phases.\n");
                    Ki_final_result = Ki_3;
                    yi_final_result(:,1) = yiv_3; % Gas
                    yi_final_result(:,2) = yiL1_3; % Liquid
                    yi_final_result(:,3) = yiL2_3; % Water
                    F_final_result(1,1) = F_3(1,1); % Gas
                    F_final_result(2,1) = Fl_3; % Liquid
                    F_final_result(3,1) = F_3(2,1); % Water
                end
        
            end
            
        else
            fprintf("\nTwo phase stability test is conducted on liquid phase composition. \n\n");
            [threephases_test_on_liquid, Ki_liquid_result, ~] = stability_test_on_two_phase_20210828(yil_2, Ki_first_stability_result, Kiv_2, p, T, Pc, Tc, omega, kij);
    
            if threephases_test_on_liquid == false 
                fprintf("\nThe two phases system is stable. \n");
                fprintf("\n Final conclusion is that the system has only two phases.\n");
        
                Ki_final_result = Kiv_2;
                yi_final_result(:,1) = yiv_2;
                yi_final_result(:,2) = yil_2;
                F_final_result(1,1) = Fv_2;
                F_final_result(2,1) = Fl_2;
        
            else
                fprintf("\nThe two phases system is not stable. \n");
                fprintf("\nSystem will split into three phases and need to perform three phase flash calculation\n");
                
                Ki_3(:,1) = 1./Kiv_2; % Two phase flash calculation
                Ki_3(:,2) = Ki_liquid_result; % Two phase stability test
                [threephase_flash, Ki_3, yiv_3, yiL2_3, yiL1_3, F_3] = vapor_liquid1_liquid2_3Peq(Ki_3, p, T, Pc, Tc, zi, omega, kij, tol, maxiter);
                Fl_3 = 1 - sum(F_3);
                
                if threephase_flash == false
                    fprintf("\n Final conclusion is that the system has only two phases.\n");
                    Ki_final_result = Kiv_2;
                    yi_final_result(:,1) = yiv_2;
                    yi_final_result(:,2) = yil_2;
                    F_final_result(1,1) = Fv_2;
                    F_final_result(2,1) = Fl_2;

                else
                    fprintf("\n Final conclusion is that the system has three phases.\n");
                    Ki_final_result = Ki_3;
                    yi_final_result(:,1) = yiv_3; % Gas
                    yi_final_result(:,2) = yiL1_3; % Liquid
                    yi_final_result(:,3) = yiL2_3; % Water
                    F_final_result(1,1) = F_3(1,1); % Gas
                    F_final_result(2,1) = Fl_3; % Liquid
                    F_final_result(3,1) = F_3(2,1); % Water
                end

            end
            
        end
        
    else
        
    fprintf("\nSystem is stable at single phase. \n");
    % For one phase.
    [~, ~] = fugacitycoef_multicomp(zi, p, T, Pc, Tc, omega, kij);
    yi_final_result = zi;
    F_final_result = 1;

    end
    
end