function non_negative_zfactor = calz_multicomp(A, B) 
%Z^3 + (B - 1)Z^2 + (A −2*B − 3*B^2)Z + (-A*B + B^2 + B^3) = 0

% Calculate the coefficients of cubic equation.
a1 = 1;
a2 = B - 1;
a3 = A - 2*B - 3*B^2;
a4 = -1*A*B + B^2 + B^3;

% Solve the cubic equation.
zroots = roots([a1 a2 a3 a4]);

% Choose the real roots.
zfactor = [];
for i = 1:length(zroots)
    if imag( zroots(i) ) == 0
        z = real(zroots(i));
        zfactor = cat(1, zfactor, z);
    end
end

non_negative_zfactor = [];
for i = 1:length(zfactor)
   if zfactor(i) >= 0
       z = zfactor(i);
       non_negative_zfactor = cat(1,non_negative_zfactor,z);    
   end
end

non_negative_zfactor = sort(non_negative_zfactor);

end