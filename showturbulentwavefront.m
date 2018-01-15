%% Final Exercise Filtering & Identification
% Niels Uitterdijk  4276892
% Jurriaan Govers   4163753

clear all
close all
clc
tic

%% Load data
% load('systemMatrices.mat')
% myfile = 'turbulenceData.mat';
% [parentdir,~,~]=fileparts(pwd);
% % load(fullfile(parentdir,myfile))
% x = phiSim{1};
load('model.mat')
%% 
phi = zeros(7,7,500);
phi(:) = yd(:);
for n = 1:5000
    surf(phi(:,:,n))
    axis([0 7 0 7 -20 20])
    pause(0.05)
end

