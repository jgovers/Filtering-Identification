function [var_eps] = AOloopMVM(G,H,C_phi0,SNR,phik);
%% System information
sig_e = sqrt(10^(-SNR/10));  
M = C_phi0*G'/(G*C_phi0*G'+sig_e^2*eye(72));
%% Initialise matrices
n = size(H,1);      % Dimension lifted wavefront
T = length(phik);   % Number of temporal phase points
eps = zeros(n,T);   % Residual wavefront
phi_h = zeros(n,T); % Residual wavefront observer
sigma = zeros(T,1); % System Variance
u = zeros(n,T);     % Control action
%% Control loop
for k = 2:T
    %% System equations
    eps(:,k) = phik(:,k) - H*u(:,k-1);
    s = awgn(G*eps(:,k),SNR);
    %% Wavefront estimator and control action
    phi_h(:,k) = M*s;
    u(:,k) = u(:,k-1) + H\phi_h(:,k);
    %% Variance
    sigma(k) = var(eps(:,k) - mean(eps(:,k)));
end
var_eps = mean(sigma);
end
