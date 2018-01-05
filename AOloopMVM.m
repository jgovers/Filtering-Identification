function [var_eps] = AOloopMVM(G,H,C_phi0,SNR,phik);
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


skIdent = awgn(G*phik,SNR);    


n = size(H,1);      % dimension lifted wavefront
% ns = size(G,1);     % dimension lifted sensor slopes
T = length(phik);   % number of temporal phase points

epsk = zeros(n,T);  % residual wavefront
eps_piston_removed = zeros(n,T); % residual wavefront with mean removed
% sk = zeros(ns,T);   % slopes measurements
% strehl = zeros(T,1);% strehl ratio
sigma = zeros(T,1);
u = zeros(n,T);


for k = 1:T-1
    epsk(:,k+1) = phik(:,k+1) - H*u(:,k);
    eps_piston_removed(:,k+1) = epsk(:,k+1) - mean(epsk(:,k+1)); 
    beta = awgn(epsk(:,k+1),0);
    u(:,k+1) = u(:,k) + H\beta;
    sigma(k+1) = var(eps_piston_removed(:,k+1));
    %     sk(:,k+1) = awgn(G*epsk(:,k+1),SNR);
    %     strehl(k+1) = exp(-sigma(k+1)^2);
end
% strehl = mean(strehl);
var_eps = mean(sigma);

end
