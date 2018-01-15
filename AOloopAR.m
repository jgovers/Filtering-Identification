function [var_eps] = AOloopAR(G,H,C_phi0,C_phi1,sig_e,phik)
%% Initialising matrices
n = size(H,1);          % State dimension 
ns = size(G,1);         % Output dimension
T = length(phik)+1;     % Number of temporal phase points
eps_h = zeros(n,T);     % Residual wavefront with stochast
eps = zeros(n,T);       % Residual wavefront with mean removed
s = zeros(ns,T);       % Output data
sigma = zeros(T,1);     % Variance of each residual wavefront measurement vector
u_k = zeros(n,T);       % Optimal control action
%% System information
[A,Cw,K] = computeKalmanAR(C_phi0,C_phi1,G,sig_e);
B = [A*H -H];
D = [G*H zeros(72,49)];
vare = -log10(sig_e)*20;
varW = 0;
%% Control loop
for k = 2:T-2
    %% System equations
    eps(:,k) = phik(:,k) - H*u_k(:,k);
    s = awgn(G*eps(:,k),vare);
    %% Observer and controller
    eps_h(:,k+1) = (A-K*G)*eps_h(:,k) + B*[u_k(:,k-1); u_k(:,k)] + K*s;
    u_k(:,k+1) = u_k(:,k) + H\eps_h(:,k+1);
    %% Variance
    sigma(k+1) = var(eps(:,k));
end
% strehl = mean(strehl);
var_eps = mean(sigma);
end
