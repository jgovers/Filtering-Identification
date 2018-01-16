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
phi = cell2mat(phiIdent);
mphi = sum(phi,2)/length(phi);
toc
% %% No Control
% fprintf('\nCalculations without control:\n')
% var_nc = zeros(size(phiIdent,2),1);
% for i = 1:size(phiIdent,2)
%     phik = phiIdent{i};
%     [var_nc(i)] = AOloop_nocontrol(phik,SNR,H,G);
%     fprintf('.')
% end
% fprintf('\n')
% toc
% %% Random walk model (3.6)
% fprintf('\nRandom Walk model:\n')
% sig_e   = sqrt(10^(-SNR/10));
% var_rw = zeros(size(phiIdent,2),1);
% for i = 1:size(phiIdent,2)
%     phik = phiIdent{i};
%     C_phi0 = cov(phik');
%     [var_rw(i)] = AOloopMVM(G,H,C_phi0,SNR,phik);
%     fprintf('.')
% end
% fprintf('\n')
% toc
% %% Kalman
% fprintf('\nKalman filter:\n')
% for i = 1:size(phiIdent,2)
%     phi = phiIdent{i};
%     C_phi0 = cov(phi');
%     C_phi1 = covariance(phi,1,mphi);
%     sig_e   = sqrt(10^(-SNR/10));
%     [var_k(i)] = AOloopAR(G,H,C_phi0,C_phi1,sig_e,phi);
%     fprintf('.')
% end
% fprintf('\n')
% toc
%% Subspace Identification
fprintf('\nSubspace Identification:\n')
Nid = 2000;
Nval = 500;
s = 25;
n = 11;
for i = 1:size(phiIdent,2)
    phi = phiIdent{i};
    [A,C,K,vaf] = N4sid(sk,Nid,Nval,s,n);
    [var_si(i)] = phiSid(G,H,A,K,C,SNR,lambda,phi);
    fprintf('.')
end
fprintf('\n')
toc
%% Plots
% figure; hold on;
% plot(var_nc); plot(var_rw); plot(var_k);
% legend('Var. Tur. Wav.','Var. Res. Wav. Random Walk model','Var. Res. Wav. Kalman filter')

