%Check phase split
function [phasesplit, Ki, Yi] = check_phase_split_on_single_phase_20210828(Test, Ki, zi, p, T, Pc, Tc, omega, kij)

syms trivial
tol = 1e-6;
maxiter = 100;
eps = 1;

[phi_zi, ~] = fugacitycoef_multicomp(zi, p, T, Pc, Tc, omega, kij);

% Initial guest of Yi
% The trial phase is gas
Yi = zi.*Ki;

for i = 1:maxiter

    % Sum of mole numbers
    S = sum(Yi);
    
    % Normalize the second-phase mole numbers to get mole fractions, yi
    % Update 20210828: Thay yi bằng Yi
    Yi = Yi./S;
    
    % Check first condition
    if max(abs(Yi - zi)) < 1e-6
        phasesplit = false;
        return;
    end
    
    if eps < tol
        break;
    end
    
    % Calculate the second-phase fugacity coefficients 
    [phi_Yi, ~] = fugacitycoef_multicomp(Yi, p, T, Pc, Tc, omega, kij);
      
    % Calculate epsilon
    eps = cal_eps(Yi, phi_Yi, zi, phi_zi);

    
    % Update new Yi
    Yi = update_Yi(zi, phi_zi, phi_Yi);
    
    Ki = Yi./zi;
    % Trivial solution check
    trivial_check = norm(log(Ki))^2;
        
    if trivial_check < 1e-4
        fprintf("First stability test (on single phase system). The trivial solution is found with the initial guest number %d at the %d iteration. \n",Test,i);
        trivial = true;
        break;
    end
    
end

if i == maxiter && S <= 1
    fprintf('The iteration in checkphasesplit() did not converge for the test number %d. One phase is assumed.\n', Test);
    phasesplit = false;
elseif S > 1
    phasesplit = true;
else
    phasesplit = false;
end

if trivial == true
    phasesplit = false;
end

%Check NaN solution
NaN_check = isnan(Ki);
if max(NaN_check) == 1
    fprintf("The solution is defined as NaN\n");
    phasesplit = false;
end

end

function eps = cal_eps(Yi, phi_Yi, zi, phi_zi)
ncomp = size(Yi,1);

eps = 0;
for i = 1:ncomp
    eps = eps + ( log(Yi(i)/zi(i)*phi_Yi(i)/phi_zi(i)) )^2;
end

end

function Yi = update_Yi(zi, phi_zi, phi_Yi)
ncomp = size(zi,1);

Yi = zeros(ncomp,1);
for i = 1:ncomp
    Yi(i) = zi(i)*phi_zi(i)/phi_Yi(i);
end

end
