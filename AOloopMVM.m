function [var_eps] = AOloopMVM(G,H,C_phi0,SNR,phik);
%% Generate measurement data
sig_e = 10^(-SNR/10);
s = awgn(G*phik,SNR);    
%% Initialise matrices
M = C_phi0*G'/(G*C_phi0*G'+sig_e^2*eye(72));
n = size(H,1);      % dimension lifted wavefront
ns = size(G,1);     % dimension lifted sensor slopes
T = length(phik);   % number of temporal phase points
eps = zeros(n,T);  % residual wavefront
eps_piston_removed = zeros(n,T); % residual wavefront with mean removed
sigma = zeros(T,1);
u = zeros(n,T);
vara = -log10(1)*10;
%% Control loop
for k = 2:T
    %% System equations
    eps(:,k) = phik(:,k) - H*u(:,k-1);
    s = awgn(G*eps(:,k),SNR);
    %% Wavefront estimator and control action
    phi_h(:,k) = M*s;
    u(:,k) = u(:,k-1) + H\phi_h(:,k);
    %% Variance
    sigma(k) = var(eps(:,k)-mean(eps(:,k)));    
end
var_eps = mean(sigma);
end
