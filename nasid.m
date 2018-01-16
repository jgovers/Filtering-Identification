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
for i=1:size(S,1)
    sv(i) = S(i,i);
end
close all
semilogy(sv,'o')
%% X estimate
Xh = sqrt(S(1:n,1:n))*V(:,1:n)';
Xh_N1 = Xh(:,1:end-1);
Xh_s1 = Xh(:,2:end);
%% System matrices
M = (Xh_N1*Xh_N1')\Xh_N1;
At = (M*Xh_s1')';
Ys1N = henk(phi,s,1,Nid-1);
Ct = (M*Ys1N')';
%% Variance matrices
W = Xh_s1 - At*Xh_N1;
V = Ys1N - Ct*Xh_N1;
Q = 1/Nid*W*W';
S = 1/Nid*W*V';
R = 1/Nid*V*V';
%% Kalman filter
[X,L,K_t] = dare(At',Ct',Q,R,S);
K = K_t';
%% Model simulation DIT IS KUT
xh = zeros(n,Nval);
yh = zeros(r,Nval);
xh(:,1) = Ct\phi(:,1);
yh(:,1) = phi(:,1);
for i = 1:Nval-1
    v = wgn(r,1,0);
    xh(:,i+1) = At*xh(:,i) + K*v;
    yh(:,i) = Ct*xh(:,i) + v;
end
%% Model validation
y = phi(:,1:Nval);
yd = y - yh;
vaf = (1-1/Nval*sum(sum(yd.*yd)))/(1/Nval*sum(sum(y.*y)))*100;
figure
plot(y(1,:))
hold on
plot(yh(1,:))
end