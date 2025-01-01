%% Computational Finance - Commodities project - A.A. 2023/2024
%% Project D - Option Calibration and Pricing (HJM) - Group 11

clear
close all
clc
format default

% Fix seed
rng(0)

addpath("questions\")
addpath("data\")
addpath("utility\")
addpath("plot\")

run("load_data.m")

%% Q3 - Model Calibration on 2026 prices

run("Q3.m");

%% Q4 - Option pricing 2026

run("Q4.m");

%% Q5 - Model Calibration on 2026 & 2028 prices

run("Q5.m");

%% Q6 - Option pricing 2026 & 2028

run("Q6.m");
