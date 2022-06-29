function [twophase_flash, Kiv, yiv, yil, Fv] = vapor_liquid_2Peq(Kiv, p, T, Pc, Tc, zi, omega, kij, tol, maxiter)

syms trivial
eps = 1;

while eps > 1e-14

    [Fv, comp] = phasefraction_multicomp(Kiv, zi, tol, maxiter);
    
    %disp(comp);

    yiv = comp(:, 1);
    yil = comp(:, 2);

    [phi_iv, ~] = fugacitycoef_multicomp_vapor_phase(yiv, p, T, Pc, Tc, omega, kij);
    [phi_il, ~] = fugacitycoef_multicomp_liquid_phase(yil, p, T, Pc, Tc, omega, kij);

    % Theo Whitson and Brule, 2000
    fli_fvi = (1./Kiv).*(phi_il./phi_iv);
    absolute_error = abs(fli_fvi - 1);
    eps = norm(absolute_error)^2;

        if eps < 1e-14
            break;
        end
    
    %Update new K
    Kiv = Kiv.*fli_fvi; 
    
    % Phuong trinh 4.51 - Whitson and Brule
    %Trivial solution check
    trivial_check = norm(log(Kiv))^2;
        if trivial_check < 1e-4
            fprintf("Trivial solution found")
            trivial = true;
            break;
        end
end

NaN_check = isnan(Kiv);
    %Check trivial solution
if trivial == true
    twophase_flash = false;
    
    %Check NaN solution
elseif max(NaN_check) == 1
    fprintf("The solution is defined as NaN\n");
    twophase_flash = false;
    
    %Check non-physical solution
elseif Fv < 0 || Fv > 1
    fprintf("The mixture is thermodynamically stable as a single phase and will not split into two phases\n");
    twophase_flash = false;
else
    twophase_flash = true;
end
end
