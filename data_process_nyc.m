clear; clc; close all;

%%%% Import Confirmed Cases
FILE ='time_series_covid19_confirmed_US.csv';
opts = detectImportOptions(FILE);
opts = setvartype(opts, opts.VariableNames, 'char');
T2 = readtable(FILE,opts);

%%%% Import Deaths
FILE ='time_series_covid19_deaths_US.csv';
opts = detectImportOptions(FILE);
opts = setvartype(opts, opts.VariableNames, 'char');
D2 = readtable(FILE,opts);

%%%% Convert to cell arrays
T = table2cell(T2);
D = table2cell(D2);

%%%% FIX: Cleanly extract the dates row before slicing columns
% This prevents indexing confusion later on.
Raw_Dates = T(1, 12:end); 

%%%% Find row indices for New York City boroughs
inds_NYC = [];
nyc_counties = {'New York', 'Kings', 'Queens', 'Bronx', 'Richmond', 'New York City'};

for i = 2:size(T,1)
    state = T{i,7};   % Column 7 is Province_State
    county = T{i,6};  % Column 6 is Admin2 (County)
    
    if isequal(state, 'New York') && ismember(county, nyc_counties)
        inds_NYC = [inds_NYC, i];
    end
end

%%%% Isolate NYC rows and strip out metadata columns immediately
% Confirmed cases has 11 metadata columns; Deaths has 12.
T_NYC = T(inds_NYC, 12:end);
D_NYC = D(inds_NYC, 13:end);

%%%% Convert from text to double-precision numbers
I_tot_NYC = zeros(size(T_NYC));
D_tot_NYC = zeros(size(D_NYC));

for i = 1:size(T_NYC,1)
    for j = 1:size(T_NYC,2)
        I_tot_NYC(i,j) = str2double(T_NYC{i,j});
        D_tot_NYC(i,j) = str2double(D_NYC{i,j});
    end
end

%%%% Aggregate the counties into a single total for NYC
I_tot_NYC = sum(I_tot_NYC, 1);
D_tot_NYC = sum(D_tot_NYC, 1);

% Calculate Daily New Cases/Deaths
I_tot_NYC_daily = I_tot_NYC(2:end) - I_tot_NYC(1:end-1);
D_tot_NYC_daily = D_tot_NYC(2:end) - D_tot_NYC(1:end-1);

% Create datamatrices for moving avg
I_tot_NYC_daily_avg = zeros(1, size(I_tot_NYC, 2));
D_tot_NYC_daily_avg = zeros(1, size(D_tot_NYC, 2));

% Loop for moving avg from 3/1/20 - 12/31/20
for d = 40:345
    I_tot_NYC_daily_avg(d) = (1/7)*sum(I_tot_NYC_daily(d-4 : d+2));
    D_tot_NYC_daily_avg(d) = (1/7)*sum(D_tot_NYC_daily(d-4 : d+2));
end

% Truncate to date range (3/1/20 - 12/31/20)
D_tot_NYC_daily_avg = max(D_tot_NYC_daily_avg(40:345), 0);
I_tot_NYC_daily_avg = max(I_tot_NYC_daily_avg(40:345), 0);
Dates = Raw_Dates(40:345);

save COVID_data_NYC_2020.mat Dates D_tot_NYC_daily_avg I_tot_NYC_daily_avg I_tot_NYC_daily D_tot_NYC_daily;
