%Có cập nhật ngày 29/8/2021

function [threephase, Ki_final_result, Yi_final_result] = stability_test_on_two_phase_20210828(zi, Ki_1p_stab, Ki_2p_flash, p, T, Pc, Tc, omega, kij)

ncomp = size(zi,1);

%Initial assumption of Ki
Ki_wilson = wilson(Pc, Tc, p, T, omega);

Ki_test(:,1) = Ki_wilson;
Ki_test(:,2) = 1./Ki_wilson;
Ki_test(:,3) = Ki_wilson.^(1/3);
Ki_test(:,4) = 1./(Ki_wilson.^(1/3));
Ki_test(:,5) = Ki_1p_stab;
Ki_test(:,6) = 1./Ki_1p_stab;
Ki_test(:,7) = Ki_2p_flash;
Ki_test(:,8) = 1./Ki_2p_flash;

% Cập nhật mới ngày 29/8/2021. Sử dụng phương trình 5.40 trong tài
% liệu (2013) REDUCED PHASE EQUILIBRIUM CALCULATIONS NEW REDUCED
% PARAMETERS, CRITICAL ANALYSIS AND FLUID CHARACTERIZATION -
% Dissertation_Gorucu_2013 folder Abbas Firoozabadi

% Hardcode - Thứ tự của CO2 luôn ở vị trí thứ 3 - sẽ cập nhật sau

% Caapjj nhật ngày 29/8/2021
i = 9;
for j = 1:ncomp
    for k = 1:ncomp
        if k == j
           Ki_test(k,i) = (0.9)/zi(k,1);
        else
           Ki_test(k,i) = (0.1/(ncomp-1))/zi(k,1);
        end
    end
    i = i + 1;
end

phase = [];
Ki = [];
Yi = [];

for i = 1:(ncomp + 8)
    test = i;
    [phase(i,1), Ki(:,i), Yi(:,i)] = check_phase_split_on_two_phase_20210828(test, Ki_test(:,i), Ki_2p_flash, zi, p, T, Pc, Tc, omega, kij);
end

Ki_result = [];
Yi_result = [];
for i = 1:(ncomp + 8)
    if phase(i,1) == true
        Ki_result = cat(2,Ki_result, Ki(:,i));
        Yi_result = cat(2,Yi_result, Yi(:,i));
    end
end

[phi_zi, ~] = fugacitycoef_multicomp(zi, p, T, Pc, Tc, omega, kij);

phi_Yi = [];
for i = 1: size(Yi_result,2)
    [phi_Yi(:,i), ~] = fugacitycoef_multicomp(Yi_result(:,i), p, T, Pc, Tc, omega, kij);
end

TPD = [];

for i = 1: size(Yi_result,2)
    TPD(i,1) = 1;
    for j = 1: ncomp
        TPD(i,1) = TPD(i,1) + Yi_result(j,i)*(log(phi_Yi(j,i)) + log(Yi_result(j,i)) - log(phi_zi(j,1)) - log(zi(j,1)) - 1);
    end
end

[~, index] = min(TPD);

Ki_final_result = Ki_result(:,index);
Yi_final_result = Yi_result(:,index);

if any(phase) == true
    threephase = true;
else
    threephase = false;
end

end
