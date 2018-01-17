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
phit = cell2mat(phiIdent);
mphi = sum(phit,2)/length(phit);
toc
%% No Control
fprintf('\nCalculations without control:\n')
var_nc = zeros(size(phiIdent,2),1);
for i = 1:size(phiIdent,2)
    phi = phiIdent{i};
    [var_nc(i)] = AOloop_nocontrol(phi,SNR,H,G);
    fprintf('.')
end
fprintf('\n')
toc
%% Random walk model (3.6)
fprintf('\nRandom Walk model:\n')
sig_e   = sqrt(10^(-SNR/10));
var_rw = zeros(size(phiIdent,2),1);
for i = 1:size(phiIdent,2)
    phi = phiIdent{i};
    C_phi0 = cov(phi');
    [var_rw(i)] = AOloopMVM(G,H,C_phi0,SNR,phi);
    fprintf('.')
end
fprintf('\n')
toc
%% Kalman filter
fprintf('\nKalman filter:\n')
% Kalman gain
C_phi0 = cov(phit');
C_phi1 = covariance(phit,1,mphi);
sig_e   = sqrt(10^(-SNR/10));
[A,Cw,K] = computeKalmanAR(C_phi0,C_phi1,G,sig_e);
% Control loop
for i = 1:size(phiIdent,2)
    phi = phiIdent{i};
    [var_k(i)] = AOloopAR(G,H,A,Cw,K,sig_e,phi);
    fprintf('.')
end
fprintf('\n')
toc
%% Subspace Identification
fprintf('\nSubspace Identification:\n')
% Settings
Nid = 3500;
Nval = 1500;
s = 17;
n = 16;
lambda = 0;
% Controller
for i = 1:size(phiIdent,2)
    phi = phiIdent{i};
    [A,C,K,vaf] = nasid(phi,Nid,Nval,s,n);
    [var_si(i)] = phiSid(G,H,A,K,C,SNR,lambda,phi);
    fprintf('.')
end
fprintf('\n')
toc
% Plots
figure; hold on;
plot(var_nc); plot(var_rw); plot(var_k);plot(var_si)
legend('Var. Tur. Wav.','Var. Res. Wav. Random Walk model','Var. Res. Wav. Kalman filter','Var. Res. Wav. Subspace Identification')

