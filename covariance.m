function [cova] = covariance(A,N,mphi)
%% Data Size
rn = length(A(:,1));
rm = length(A(1,:));
%% Calculation Covariance Matrix
cova = zeros(min(rn,rm),min(rn,rm));
for n = 2:rm
    Ci = (A(:,n)-mphi)*(A(:,n-N)-mphi)';
    cova = cova + Ci;
end
cova = cova/rm;
end
