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
% function [ sigma ] = AOloop_nocontrol(phik,SNR,H,G)



[sigma] = AOloop_nocontrol(phiSim{1,1},SNR,H,G);

%% Test test