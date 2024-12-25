clear all
close all
clc

%% Read Prices
path_map = strcat(pwd, '/data/');
price_filename = 'prices.xlsx';

table_prices = readtable(strcat(path_map, price_filename));
%% Transform prices from table to timetable
dt = table_prices(:, 1).Variables;
values = table_prices(:, 2:end).Variables;
assetNames = table_prices.Properties.VariableNames(2:end);
nAssets = length(assetNames);

myPrice_dt = array2timetable(values, 'RowTimes', dt, 'VariableNames', assetNames);
%% Selection of a subset of Dates
start_dt = datetime('01/01/2024', 'InputFormat', 'dd/MM/yyyy');
end_dt = dt(end);

rng = timerange(start_dt, end_dt, 'closed');
subsample = myPrice_dt(rng, :);

prices_val = subsample.Variables;
dates_ = subsample.Time;

%% Calculate returns
Ret = prices_val(2:end, :)./prices_val(1:end-1, :);
LogRet = log(Ret);
ExpRet = mean(LogRet);
CovMatrix = cov(LogRet);

%% Load portfolios and put them in a 'portfolios' object
load('portfolios.mat')

allObjects = who;
portfolios_names = allObjects(contains(allObjects, 'portfolio'));

portfolios = cell(1, length(portfolios_names));
for i = 1:length(portfolios_names)
     portfolios{i} = eval(portfolios_names{i});
end

%% OOS evaluation
oos_evals = cell(1, length(portfolios)+1);
for i = 1:length(portfolios)
    oos_evals{i} = {portfolios{i}.evaluate(Ret, 'saveResults', false, 'printTable', true)};
end


PortfolioViz.plotEvolution(portfolios, dates_, Ret);
title("Out of sample performance of all the portfolios");


