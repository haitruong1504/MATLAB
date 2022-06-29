function g = cal_gibbsenergy(yi, fugcoef, p)

ncomp = size(yi,1);
g = 0;

for i = 1:ncomp
    if yi(i) ~= 0
        g = g + yi(i)*log(yi(i)*fugcoef(i)*p);
    end
end
