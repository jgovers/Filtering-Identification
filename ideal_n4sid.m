function [Nid,s,n,minvaf] = ideal_n4sid(phit,Nval,x)
Nid = 4000;
for j=1:x
    parfor k=1:x
        n = 8+k;
        s = n+j;
        [A,C,K,vaf(j,k)] = nasid(phit,Nid,Nval,s,n);
    end
fprintf('.')
end
[minM idx] = max(vaf(:));
[xc yc] = ind2sub(size(vaf),idx);
n = 8+yc;
s = n+xc;
minvaf = minM;
end