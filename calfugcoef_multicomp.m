function fugcoef = calfugcoef_multicomp(zfactor, A, B, Am, Bm, Am2)

ncomp = size(A,1);
fugcoef = zeros(ncomp,1);

c0 = 2*sqrt(2);
c1 = 1 + sqrt(2);
c2 = 1 - sqrt(2);

for i = 1:ncomp
    
    if zfactor < Bm
        error('Z-factor must be larger than Bmix.\n');
    end
    
    fugcoef(i) = (B(i)/Bm)*(zfactor - 1) - log(zfactor - Bm) ...
        - (Am/(c0*Bm))*(2*Am2(i)/Am - B(i)/Bm) ...
        *log((zfactor + c1*Bm)/(zfactor + c2*Bm));
    
    fugcoef(i) = exp(fugcoef(i));
end

end