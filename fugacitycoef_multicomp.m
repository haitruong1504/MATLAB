%% CALCULATE THE FUGACITY COEFFICIENT AND Z-FACTOR OF MULTI-COMPONENT SYSTEMS
% -------------------------------------------------------------------------
% The Definition of Variables.
% comp    : composition
% press   : pressure
% temp    : temperature
% pressc  : critical pressure
% tempc   : critical temperature
% acentric: acentric factor
% BIP     : binary interaction parameter
% fugcoef : fugacity coefficient
% zfactor : compressibility factor
% -------------------------------------------------------------------------
% In this function, an appropriate z-factor is automatically chosen
% according to gibbs free energy if multiple roots are found.
function [fugcoef, zfactor] = fugacitycoef_multicomp(yi, p, T, Pc, Tc, omega, kij)

[Ai, Bi] = calAB_multicomp(omega, p, T, Pc, Tc);

[Am, Bm, Am2] = calABmixture_multicomp(Ai, Bi, yi, kij);

zfactor = calz_multicomp(Am, Bm);

if (size(zfactor,1) > 1)
    zfactor = choosezfactor(zfactor, p, yi, Ai, Bi, Am, Bm, Am2);
end

fugcoef = calfugcoef_multicomp(zfactor, Ai, Bi, Am, Bm, Am2);

end

%% SEARTCH AND RETURN AN APPROPRIATE Z-FACTOR
% Calculate dimensionless excess gibbs free energy, and return the z
% factor which minimizes the gibbs free energy.
function minzfactor = choosezfactor(zfactor, p, yi, Ai, Bi, Am, Bm, Am2)

gibbsenergy = [];

for i = 1:size(zfactor,1)
    fugcoef = calfugcoef_multicomp(zfactor(i), Ai, Bi, Am, Bm, Am2);
    g = calcgibbsenergy(yi, fugcoef, p);
    gibbsenergy = cat(1,gibbsenergy,g);
end

[~, index] = sort(gibbsenergy);
minzfactor = zfactor(index(1));

end

function g = calcgibbsenergy(yi, fugcoef, p)

ncomp = size(yi,1);
g = 0;

for i = 1:ncomp
    if yi(i) ~= 0
        g = g + yi(i)*log(yi(i)*fugcoef(i)*p);
    end
end

end