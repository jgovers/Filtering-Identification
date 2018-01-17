function [At,Ct,K,vaf] = nasid(phi,Nid,Nval,s,n)
%% Data size
r = size(phi,1);
%% Henkel matrices
YssN = henk(phi,s,s,Nid);
Y0sN = henk(phi,1,s,Nid);
%% QR factorization
[Qt,Rt] = qr([Y0sN;YssN]');
Q = Qt'; R = Rt';
R11 = R(1:r*s,1:r*s);
R21 = R(r*s+1:2*r*s,1:r*s);
%% SVD factorization
[U,S,V] = svd(R21/R11*Y0sN);
% for ii = 1:size(S,1)
%     Sd(ii) = S(ii,ii);
% end
% figure
% semilogy(Sd,'o')
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
Qc = 1/Nid*W*W';
Sc = 1/Nid*W*V';
Rc = 1/Nid*V*V';
%% Kalman filter
[X,L,K_t] = dare(At',Ct',Qc,Rc,Sc);
K = K_t';
%% Model simulation
xh = zeros(n,Nval);
yh = zeros(r,Nval);
xh(:,1) = Ct\phi(:,1);
yh(:,1) = phi(:,1);
AKC = At-K*Ct;
for i = 1:Nval-1
    xh(:,i+1) = AKC*xh(:,i) + K*phi(:,i);
    yh(:,i) = Ct*xh(:,i);
end
%% Model validation
vaf = (1-(var(phi(:,1:Nval) - yh)/var(phi(:,1:Nval))))*100;
end