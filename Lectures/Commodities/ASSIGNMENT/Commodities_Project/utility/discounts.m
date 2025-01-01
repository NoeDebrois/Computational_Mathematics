function [discounts, rates] = discounts(dat,dis,t0,dates)

%--------------------------------------------------------------------------

% Compute discounts and interest rates for a given set of dates

% INPUT
% dat
% dis
% t0
% dates >> datetime array of dates

% OUTPUT
% discounts
% rates

%--------------------------------------------------------------------------

discounts = interp1(dat, dis, exceltime(dates));
[~,ny] = ndays(t0, dates);
rates = -log(discounts) ./ ny; 

end

