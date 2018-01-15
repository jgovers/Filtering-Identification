function [A,Cw,K] = computeKalmanAR(C_phi0,C_phi1,G,sig_e)
A = C_phi1/C_phi0;
Cw = triu(C_phi0 - A*C_phi0*A');
Cw = Cw + Cw';
[X,L,K_t] = dare(A',G',Cw,sig_e^4*eye(72));
K = K_t';
end