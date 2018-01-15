function [var_eps] = AOloopAR(G,H,C_phi0,C_phi1,sig_e,phik)
%% Initialising matrices
n = size(H,1);          % State dimension 
ns = size(G,1);         % Output dimension
T = length(phik);     % Number of temporal phase points
eps_h = zeros(n,T);     % Residual wavefront with stochast
eps = zeros(n,T);       % Residual wavefront with mean removed
sigma = zeros(T,1);     % Variance of each residual wavefront measurement vector
u_k = zeros(n,T);       % Optimal control action
%% System information
[A,Cw,K] = computeKalmanAR(C_phi0,C_phi1,G,sig_e);
vare = -log10(sig_e)*20;
%% Control loop
for k = 2:T
    %% System equations
    eps(:,k) = phik(:,k) - H*u_k(:,k-1);
    s = awgn(G*eps(:,k),vare);
    %% Observer and controller
    u_k(:,k) = H\A*(H*u_k(:,k-1) + eps_h(:,k));
    eps_h(:,k+1) = (A-K*G)*eps_h(:,k) + A*H*u_k(:,k-1) - H*u_k(:,k) + K*s;
end
mean_eps = sum(eps(2:end),2)/T;
var_eps = sum((eps(2:end) - mean_eps).^2,2)/T;
var_eps = sqrt(sum(var_eps)/T);
end
