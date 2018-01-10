function [var_eps] = AOloopMVM(G,H,C_phi0,sig_e,phik);
%% Generate measurement data
SNR = -log10(sig_e)*10;
s = awgn(G*phik,SNR);    
%% Initialise matrices
n = size(H,1);      % dimension lifted wavefront
ns = size(G,1);     % dimension lifted sensor slopes
T = length(phik);   % number of temporal phase points
eps = zeros(n,T);  % residual wavefront
eps_piston_removed = zeros(n,T); % residual wavefront with mean removed
sigma = zeros(T,1);
u = zeros(n,T);
vara = -log10(1)*10;
%% Control loop
for k = 1:T-1
    %% Wavefront estimator and control action
    phi_h(:,k) = C_phi0*G'*(G*C_phi0*G'+sig_e^2)^(-1)*s(:,k);
    u(:,k+1) = phi_h(:,k);
    %% System response
    eps(:,k+1) = phik(:,k+1) - H*u(:,k);
    sigma(k+1) = var(eps(:,k+1)-mean(eps(:,k+1)));    
end
var_eps = mean(sigma);
end
