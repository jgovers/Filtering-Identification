%% Final Exercise Filtering & Identification
% Niels Uitterdijk  4276892
% Jurriaan Govers   4163753

clear all
close all
clc
tic

%% Load data
load('systemMatrices.mat')
 myfile = 'turbulenceData.mat';
 [parentdir,~,~]=fileparts(pwd);
 load(fullfile(parentdir,myfile))

 
%% 3.6
toc
sig_e   = sqrt(10^(-SNR/10));

var_eps = zeros(size(phiIdent,2),1);

for i = 1:size(phiIdent,2)
    phik = phiIdent{i};
    C_phi0 = cov(phik');
    [var_eps(i)] = AOloopMVM(G,H,C_phi0,SNR,phik);
end

%% No Control
toc
sigma = zeros(size(skIdent,2),1);

for i = 1:size(skIdent,2)
    phik = skIdent{i};
    [sigma(i)] = AOloop_nocontrol(phik,SNR,H,G);
end

%% Kalman
phik = phiIdent{1};
C_phi0 = cov(phik');
C_phi1 = covariance(phik,1);
sig_e   = sqrt(10^(-SNR/10));
[A,Cw,K] = computeKalmanAR(C_phi0,C_phi1,G,sig_e);



%% Plots
toc
figure; hold on;
plot(sigma); plot(var_eps);
legend('Sigma no control','var control')

