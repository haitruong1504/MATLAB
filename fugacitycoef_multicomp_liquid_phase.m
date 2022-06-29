%% CALCULATE THE FUGACITY COEFFICIENT AND Z-FACTOR OF MULTI-COMPONENT SYSTEMS
% -------------------------------------------------------------------------
% The Definition of Variables.
% yi      : composition
% p       : pressure
% T       : temperature
% Pc      : critical pressure
% Tc      : critical temperature
% omega: acentric factor
% kij     : binary interaction parameter
% fugcoef : fugacity coefficient
% zfactor : compressibility factor
% -------------------------------------------------------------------------
% In this function, the minimum z-factor is automatically chosen.
function [fugcoef, zfactor] = fugacitycoef_multicomp_liquid_phase(yi, p, T, Pc, Tc, omega, kij)

[Ai, Bi] = calAB_multicomp(omega, p, T, Pc, Tc);

[Am, Bm, Am2] = calABmixture_multicomp(Ai, Bi, yi, kij);

zfactor = calz_multicomp(Am, Bm);

if (size(zfactor,1) > 1)
    zfactor = min(zfactor);
end

fugcoef = calfugcoef_multicomp(zfactor, Ai, Bi, Am, Bm, Am2);

end
