function [ sigma ] = AOloop_nocontrol(phik,SNR,H,G)
%% Initialising matrices
n = size(H,1);      % dimension lifted wavefront
ns = size(G,1);     % dimension lifted sensor slopes
T = length(phik);   % number of temporal phase points
epsk = zeros(n,T);  % residual wavefront
sigma = zeros(T,1);
for k = 1:T-1
    epsk(:,k+1) = phik(:,k+1);
    sigma(k+1) = var(epsk(:,k+1)-mean(epsk(:,k+1));
end
sigma = mean(sigma);
end



