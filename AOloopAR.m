function [var_eps] = AOloopAR(G,H,C_phi0,C_phi1,sig_e,A,Cw,K,phik)
%% Wat te doen
% 0. Find out how to add white noise with covariance MATRIX
% 1. Generate s(k) from phi(k) by s(k) = G phi(k) + e(k)~(0,sigma_e^2).
% epsh(k+1) = C_phi\C_phi1*eps(k) + [AH -H][u(k-1) u(k)]' + K(y(k) - C eps(k))

% eps(k+1) = C_phi\C_phi1*eps(k) + [AH -H][u(k-1) u(k)]' + w(k)
% s(k) = G eps(k) + [GH 0][u(k-1) u(k)]' + e(k)
%% Initialising matrices
n = size(H,1);          % State dimension 
ns = size(G,1);         % Output dimension
T = length(phik)+1;     % Number of temporal phase points
eps_h = zeros(n,T);     % Residual wavefront with stochast
eps = zeros(n,T);       % Residual wavefront with mean removed
sk = zeros(ns,T);       % Output data
sigma = zeros(T,1);     % Variance of each residual wavefront measurement vector
u_k = zeros(n,T);       % Optimal control action
%% System information
B = [A*H -H];
D = [G*H zeros(72,49)];
vare = -log10(sig_e)*10;
varW = 0;
%% Generation of y data
sk = awgn(G*phik,vare);
%% Eerst epshat, dan gebruiken voor u. vervolgens var berekenen van var van epshat.
for k = 2:T-2
    eps_h(:,k+1) = A*eps_h(:,k) + B*[u_k(:,k-1); u_k(:,k)] + K*(sk(:,k) - G*eps_h(:,k));
    u_k(:,k+1) = u_k(:,k) + H\eps_h(:,k+1);
    eps(:,k) = phik(:,k) - H*u_k(:,k);
    sigma(k+1) = var(eps(:,k));
end
% strehl = mean(strehl);
var_eps = mean(sigma);
end
