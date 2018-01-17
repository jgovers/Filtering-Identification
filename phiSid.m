function [var_si] = phiSid(G,H,A,K,C,SNR,lambda,phi)
%% Data size
[s,l] = size(phi);
n = size(A,1);
%% Initialise matrices
eps = zeros(s,l);   % Residual wavefront
phi_h = zeros(s,l); % Residual wavefront observer
sigma = zeros(l,1); % System Variance
u = zeros(s,l);     % Control action
x = zeros(n,l);     % Wavefront states
%% Control loop
M = (H'*H + lambda*eye(s))\H';
AKC = A-K*C;
for i=2:l
    %% System equation
    eps(:,i) = phi(:,i) - H*u(:,i-1);
    %% Wavefront observer
    x(:,i+1) = AKC*x(:,i) + K*phi(:,i);
    phi_h(:,i+1) = C*x(:,i+1);
    %% Control action
    u(:,i) = M*phi_h(:,i+1);
    %% Variance
    sigma(i) = var(eps(:,i) - mean(eps(:,i)));
end
var_si = mean(sigma);
end