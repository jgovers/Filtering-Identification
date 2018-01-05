function [A,Cw,K] = computeKalmanAR(C_phi0,C_phi1,G,sig_e)

A = C_phi0\C_phi1;
Cw = C_phi0 - A*C_phi0*A';
[X,L,K_t] = dare(A',eye(49),0,Cw);
K = K_t';

end