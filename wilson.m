function K = wilson(Pc, Tc, p, T, omega)

K  = zeros(length(omega),1);% equilibrium constant

% estimate K by using Wilson equation.
for i = 1:length(omega)
    
    K(i) = (Pc(i)/p)*(exp(5.373*(1 + omega(i))*(1 - Tc(i)/T)));
    
end

end

