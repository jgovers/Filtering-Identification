function [var_eps] = AOloopAR(G,H,C_phi,sig_e,A,Cw,phi)
% Title of function
% IN
% G         : measurement matrix 
% H         : influence matrix mapping the wavefront on the mirror
% C_phi0    : variance matrix
% sig_e     : variance
% phik      : incoming turbulence wavefront

% SNR  : Signal to Noise Ratio for the sensor


% OUT
% var_eps   : mean variance of the residual wavefront


n = size(H,1);      % dimension lifted wavefront
% ns = size(G,1);     % dimension lifted sensor slopes
T = length(phik)+1;   % number of temporal phase points

epskhat = zeros(n,T);  % residual wavefront
eps_piston_removed = zeros(n,T); % residual wavefront with mean removed
% sk = zeros(ns,T);   % slopes measurements
% strehl = zeros(T,1);% strehl ratio
sigma = zeros(T,1);
u = zeros(n,T);


for k = 2:T-1
    epshatk(:,k+1) = A*epshatk(:,k) + [A*H,-H]*[u(:,k-1),u(:,k)] + 
    eps_piston_removed(:,k+1) = epskhat(:,k+1) - mean(epskhat(:,k+1)); 
    beta = awgn(epskhat(:,k+1),0);
    u(:,k+1) = u(:,k) + H\beta;
    sigma(k+1) = var(eps_piston_removed(:,k+1));
    %     sk(:,k+1) = awgn(G*epsk(:,k+1),SNR);
    %     strehl(k+1) = exp(-sigma(k+1)^2);
end
% strehl = mean(strehl);
var_eps = mean(sigma);


end
