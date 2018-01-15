function [HssN] = henk(data,i,s,N)
HssN = zeros(s,N);
for n = 1:s
    for m = 1:N
        HssN(n,m) = data(n+m);
    end
end
end