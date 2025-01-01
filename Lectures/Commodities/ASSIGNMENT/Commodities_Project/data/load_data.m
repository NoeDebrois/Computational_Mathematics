clear

path_map = strcat(pwd, '/data/');
price_filename = 'DATA_FREEX.xlsx';
path = strcat(path_map, price_filename);

data.fwd = readtable(path, 'Sheet', 'Prices');

% Read the 2026 data
data2026 = readtable(path, 'Sheet', '2026 surface');
volsurface_2026.strikes = table2array(data2026(1, 2:end))';
volsurface_2026.tenors = table2array(data2026(2:end, 1));
volsurface_2026.vols = table2array(data2026(2:end, 2:end));
data.volsurface_2026 = volsurface_2026;

% Read the 2028 data
data2028 = readtable(path, 'Sheet', '2028 surface');
volsurface_2028.strikes = table2array(data2028(1, 2:end))';
volsurface_2028.tenors = table2array(data2028(2:end, 1));
volsurface_2028.vols = table2array(data2028(2:end, 2:end));
data.volsurface_2028 = volsurface_2028; 

% Discount data
data.discounts = table2array(readtable(path, 'Sheet', 'Discounts'));

% Initial time
data.t0 = datetime('04-Nov-2024',"InputFormat","dd-MMM-uuuu");

clear data2026 data2028 fwd path path_map price_filename 
clear volsurface_2026 volsurface_2028 

disp('> Data loaded')