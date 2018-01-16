function [HssN] = henk(data,o,s,N)
[n,m] = size(data);
HssN = zeros(s*n,N);
for i = 1:s
    for j = 1:N
        HssN(n*i-n+1:n*i,j) = data(:,i+j+o-2);
    end
end
end