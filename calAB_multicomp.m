function [Ai, Bi] = calAB_multicomp(omega, p, T, Pc, Tc)
%Constant:
R = 10.7315; %psia*ft3/((lbm.mol)*R)
omega_a_0 = 0.45724;
omega_b_0 = 0.07780;

Tr = T./Tc;

%Calculate intermediate variable for estimating Z
mi = zeros (length(omega),1);
alphai = zeros (length(omega),1);
ai = zeros (length(omega),1);
bi = zeros (length(omega),1);
Ai = zeros (length(omega),1);
Bi = zeros (length(omega),1);
for i = 1: length(omega)
    if omega(i) > 0.49
        mi(i) = 0.3796 + 1.485*omega(i) - 0.1644*(omega(i))^2 + 0.01667*(omega(i))^3;
    else 
        mi(i) = 0.37464 + 1.54226*omega(i) - 0.26992*(omega(i))^2;
    end
alphai(i) = (1 + mi(i)*(1- sqrt(Tr(i))))^2;

ai(i) = omega_a_0*(((R^2)*(Tc(i)^2))/Pc(i))*alphai(i);
bi(i) = omega_b_0*((R*Tc(i))/Pc(i));

Ai(i) = (ai(i)*p)/((R*T)^2);
Bi(i) = (bi(i)*p)/(R*T);
end

end
