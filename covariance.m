function [cova] = covariance(A,N)
%% Size of A
rn = length(A(:,1));
rm = length(A(1,:));
%% Check if orientation of A is alright
[r,l] = min([rn rm]);
if l == 2
    A = A';
end
%% Calculation of variance matrix
cova = zeros(r,r);
for n=1:r
    for m = 1:r
        if n == m
            cova_t = cov(A(n,1:rm-N)-mean(A(n,1:rm-N)),A(n,N+1:rm)-mean(A(n,N+1:rm)));
            cova(n,m) = cova_t(2);
        else
            cova_t = cov(A(n,1:rm-N)-mean(A(n,1:rm-N)),A(m,N+1:rm)-mean(A(m,N+1:rm)));            
            cova(n,m) = cova_t(2);
        end
    end
end
if l == 2
    cova = cova';
end
fprintf('This variance matrix is generated specificly for Filtering & Identification!\n\n')
end
