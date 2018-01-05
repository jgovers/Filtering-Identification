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
C_phi0  = 0;
sig_e   = 0;
% phi     = 0;
var_eps = zeros(size(phiSim,2),1);

for i = 1:size(phiSim,2)
    phik = phiSim{i};
    [var_eps(i)] = AOloopMVM(G,H,C_phi0,sig_e,phik);
end

%% No Control
toc
sigma = zeros(size(phiSim,2),1);

for i = 1:size(phiSim,2)
    phik = phiSim{i};
    [sigma(i)] = AOloop_nocontrol(phik,SNR,H,G);
end

%% Kalman
toc


%% Plots
toc
figure; hold on;
plot(sigma); plot(var_eps);
legend('Sigma no control','var control')

