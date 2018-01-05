%% Final Exercise Filtering & Identification
% Niels Uitterdijk  4276892
% Jurriaan Govers   4163753

clear all
close all
clc

%% Load data
load('systemMatrices.mat')
 myfile = 'turbulenceData.mat';
 [parentdir,~,~]=fileparts(pwd);
 load(fullfile(parentdir,myfile))
 
%% 3.6
C_phi0  = 0;
sig_e   = 0;
phi     = 0;

[var_eps] = AOloopMVM(G,H,C_phi0,sig_e,phi);

%% No Control
sigma = zeros(size(phiSim));
for i = 1:size(phiSim,2)
    phik = phiSim{i};
    [sigma(i)] = AOloop_nocontrol(phik,SNR,H,G);
end

plot(sigma)
