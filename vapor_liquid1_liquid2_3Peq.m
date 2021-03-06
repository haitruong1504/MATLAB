function [threephase_flash, Ki, yiv, yiL2, yiL1, F] = vapor_liquid1_liquid2_3Peq(Ki, p, T, Pc, Tc, zi, omega, kij, tol, maxiter)

syms trivial 
syms negative_flash
eps1 = 1;
eps2 = 1;

while eps1 > 1e-13 && eps2 > 1e-13

    [F, comp] = phasefraction_multicomp(Ki, zi, tol, maxiter);

    yiv = comp(:,1);
    yiL2 = comp(:,2);
    yiL1 = comp(:,3);

    [phi_iv, ~] = fugacitycoef_multicomp_vapor_phase(yiv, p, T, Pc, Tc, omega, kij);
    
    [phi_iL2, ~] = fugacitycoef_multicomp_liquid_phase(yiL2, p, T, Pc, Tc, omega, kij);
    [phi_iL1, ~] = fugacitycoef_multicomp_liquid_phase(yiL1, p, T, Pc, Tc, omega, kij);

    fL1i_fvi = (1./Ki(:,1)).*(phi_iL1./phi_iv);
    absolute_error_1 = abs(fL1i_fvi - 1);
    eps1 = norm(absolute_error_1)^2;
    
    fL1i_fL2i = (1./Ki(:,2)).*(phi_iL1./phi_iL2);
    absolute_error_2 = abs(fL1i_fL2i - 1);
    eps2 = norm(absolute_error_2)^2;

        if eps1 < 1e-13 && eps2 < 1e-13
            break;
        end

    Ki(:,1) = Ki(:,1).*fL1i_fvi;
    Ki(:,2) = Ki(:,2).*fL1i_fL2i;
    
    trivial_check_1 = norm(log(Ki(:,1)))^2;
    trivial_check_2 = norm(log(Ki(:,2)))^2;
        if trivial_check_1 < 1e-4 && trivial_check_2 < 1e-4
            fprintf("Trivial solution found")
            trivial = true;
            break;
        end
end

for i = 1:length(F)
    if F(i) > 1 || F(i) < 0
        fprintf("The mixture is thermodynamically stable as a two phase and will not split into three phases.\n");
        negative_flash = true;
    end
end

NaN_check = isnan(Ki);
    %Check trivial solution
if trivial == true
    threephase_flash = false;
    
    %Check NaN solution
elseif max(NaN_check) == 1
    fprintf("The solution is defined as NaN\n");
    threephase_flash = false;
    
    %Check non-physical solution
elseif negative_flash == true
    
    threephase_flash = false;
else
    threephase_flash = true;
end


