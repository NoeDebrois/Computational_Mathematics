function [n, ny] = ndays(t0, t1)

%--------------------------------------------------------------------------

% Retrive number of days between 2 datetimes

% INPUT
% t0
% t1

% OUTPUT
% n >> Number of days
% ny >> Number of days as a fraction of a year

%--------------------------------------------------------------------------

n = caldays(between(t0, t1, "days"));
ny = n/365;

end

