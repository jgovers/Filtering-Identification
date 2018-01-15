function [var_eps] = AOloopAR(G,H,C_phi0,C_phi1,sig_e,phik)
%% Initialising Matrices
n = size(H,1);          % State dimension 
ns = size(G,1);         % Output dimension
T = length(phik);     % Number of temporal phase points
eps_h = zeros(n,T);     % Residual wavefront with stochast
eps = zeros(n,T);       % Residual wavefront with mean removed
sigma = zeros(T,1);     % Variance of each residual wavefront measurement vector
u_k = zeros(n,T);       % Optimal control action
%% System Information
vare = -log10(sig_e)*20;
[A,Cw,K] = computeKalmanAR(C_phi0,C_phi1,G,sig_e);
%% Control Loop
for k = 2:T
    %% System Equations
    eps(:,k) = phik(:,k) - H*u_k(:,k-1);
    s = awgn(G*eps(:,k),vare);
    %% Observer and Controller
    u_k(:,k) = H\((A-K*G)*eps_h(:,k)+K*s + A*H*u_k(:,k-1));
    eps_h(:,k+1) = (A-K*G)*eps_h(:,k) + A*H*u_k(:,k-1) - H*u_k(:,k) + K*s;
    %% Variance
    sigma(k) = var(eps(:,k) - mean(eps(:,k)));
end
var_eps = mean(sigma);
end
