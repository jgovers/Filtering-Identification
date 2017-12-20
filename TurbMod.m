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
%% 3.1 Reconstruct the waveform
