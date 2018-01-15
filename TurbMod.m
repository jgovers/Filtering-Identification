%% Final Exercise Filtering & Identification
% Niels Uitterdijk  4276892
% Jurriaan Govers   4163753
clear all
close all
clc
tic
%% Load data
fprintf('Loading data:\n')
load('systemMatrices.mat')
myfile = 'turbulenceData.mat';
[parentdir,~,~]=fileparts(pwd);
load(fullfile(parentdir,myfile))
toc
% %% No Control
% fprintf('\nCalculations without control:\n')
% var_nc = zeros(size(phiIdent,2),1);
% for i = 1:size(phiIdent,2)
%     phik = phiIdent{i};
%     [var_nc(i)] = AOloop_nocontrol(phik,SNR,H,G);
% end
% toc
% %% Random walk model (3.6)
% fprintf('\nRandom Walk model:\n')
% sig_e   = sqrt(10^(-SNR/10));
% var_rw = zeros(size(phiIdent,2),1);
% for i = 1:size(phiIdent,2)
%     phik = phiIdent{i};
%     C_phi0 = cov(phik');
%     [var_rw(i)] = AOloopMVM(G,H,C_phi0,SNR,phik);
% end
% toc
%% Kalman
fprintf('\nKalman filter:\n')
for i = 1:size(phiIdent,2)
    phik = phiIdent{i};
    C_phi0 = cov(phik');
    C_phi1 = covariance(phik,1);
    sig_e   = sqrt(10^(-SNR/10));
    [var_k(i)] = AOloopAR(G,H,C_phi0,C_phi1,sig_e,phik);
end
toc
%% Plots
figure; hold on;
plot(var_nc); plot(var_rw);
legend('var no control','var control')

