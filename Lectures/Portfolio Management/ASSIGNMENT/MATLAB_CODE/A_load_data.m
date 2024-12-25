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
start_dt = datetime('01/01/2023', 'InputFormat', 'dd/MM/yyyy');
end_dt = datetime('31/12/2023', 'InputFormat', 'dd/MM/yyyy');

rng = timerange(start_dt, end_dt, 'closed');
subsample = myPrice_dt(rng, :);

prices_val = subsample.Variables;
dates_ = subsample.Time;

%% Define Sector Assignments for Assets
sectors = struct();
sectors.cyclical = {'ConsumerDiscretionary', 'Financials', 'Materials', 'RealEstate', 'Industrials'};
sectors.defensive = {'ConsumerStaples', 'Utilities', 'HealthCare'};
sectors.sensible = {'Energy', 'InformationTechnology', 'CommunicationServices'};

% Get indices for sector assignments
cyclical_indices = double(ismember(assetNames, sectors.cyclical));
defensive_indices = double(ismember(assetNames, sectors.defensive));
sensible_indices = double(ismember(assetNames, sectors.sensible));

%% Calculate returns
Ret = prices_val(2:end, :)./prices_val(1:end-1, :);
LogRet = log(Ret);
ExpRet = mean(LogRet);
CovMatrix = cov(LogRet);

%% Read Capitalizations
cap_filename = 'capitalizations.xlsx';

table_caps = readtable(strcat(path_map, cap_filename));
caps = table_caps(:, 2:end).Variables;

disp("A_load_data.m > Data from prices.xslx and capitalization.xlsx loaded.");
disp(" ");