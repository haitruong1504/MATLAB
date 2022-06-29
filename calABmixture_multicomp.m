function [Am, Bm, Am2] = calABmixture_multicomp(Ai, Bi, yi, kij)

Aij = zeros(length(yi), length(yi));
for i = 1:length(yi)
    for j = 1:length(yi)
        Aij(i,j) = sqrt(Ai(i)*Ai(j))*(1 - kij(i,j));
    end
end

Am = 0;
for i = 1:length(yi)
    for j = 1:length(yi)
    Am = Am + yi(i)*yi(j)*Aij(i,j);
    end
end

Bm = 0;
for i = 1:length(yi) 
    Bm = Bm + yi(i)*Bi(i);
end

Am2 = zeros(length(yi),1);
for i = 1:length(yi)
    for j = 1:length(yi)
        Am2(i) = Am2(i) + yi(j)*Aij(i,j);
    end
end

end
