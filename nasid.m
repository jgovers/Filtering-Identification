function [At,Ct,K,vaf] = nasid(phi,Nid,Nval,s,n)
%% Data properties
[r,m] = size(phi);
%% Henkel matrices
YssN = henk(phi,s,s,Nid);
Y0sN = henk(phi,0,s,Nid);
%% QR factorization
[Qt,Rt] = qr([Y0sN;YssN]');
Q = Qt'; R = Rt';
R22 = R(1:r*s,1:r*s);
R32 = R(r*s+1:2*r*s,1:r*s);
%% SVD factorization
[U,S,V] = svd(R32/R22*Y0sN);
for i=1:s%size(S,1)
    sv(i) = S(i,i);
end
semilogy(sv,'o')
%% X estimate
Xh = chol(S(1:n,1:n))*V(:,1:n)';
Xh_N1 = Xh(:,1:end-1);
Xh_s1 = Xh(:,2:end);
%% System matrices DIT IS KUT
M = (Xh_N1*Xh_N1')\Xh_N1;
At = (M*Xh_s1')';
Ys1N = henk(phi,s,1,Nid-1);
Ct = (M*Ys1N')';
%% Variance matrices
Wh = Xh_s1 - At*Xh_N1;
Vh = Ys1N - Ct*Xh_N1;
Qh = 1/Nid*(Wh*Wh');
Sh = 1/Nid*(Wh*Vh');
Rh = 1/Nid*(Vh*Vh');
%% Kalman filter
[X,L,K_t] = dare(At',Ct',Qh,Rh,Sh);
K = K_t';
%% Model simulation DIT IS KUT
xh = zeros(n,Nval);
yh = zeros(r,Nval);
xh(:,1) = Ct\phi(:,1);
yh(:,1) = phi(:,1);
for i = 1:Nval-1
    v = randn(49,1);
    xh(:,i+1) = At*xh(:,i) + K*v;
    yh(:,i) = Ct*xh(:,i) + v;
end
%% Model validation
y = phi(:,1:Nval);
yd = y - yh;
vaf = max(0,(1-1/Nval*sum(sum(yd.*yd))/(1/Nval*sum(sum(y.*y))))*100);
end