cd("C:\Users\Italo\OneDrive\Documents\PhD\PhD\First")
addpath('C:\Users\Italo\OneDrive\Documents\PhD\PhD\First\ddisk\matlab');
addpath('C:\Users\Italo\OneDrive\Documents\PhD\PhD\First\ddisk\m_utilities');
addpath('C:\Users\Italo\OneDrive\Documents\PhD\PhD');
addpath('C:\Users\Italo\OneDrive\Documents\PhD\PhD\First\FFR_code_final\');
addpath('C:\Users\Italo\OneDrive\Documents\PhD\PhD\First\FFR_code_final\functions\ML');
addpath('C:\Users\Italo\OneDrive\Documents\PhD\PhD\First\FFR_code_final\functions\Benchmarks\');
addpath('C:\Users\Italo\OneDrive\Documents\PhD\PhD\First\FFR_code_final\functions\PC');

%Retrieveing yield data
yield_data=readtable("C:\Users\Italo\OneDrive\Documents\PhD\PhD\USD Zero Coupon Yields - Nice names.xlsx");

yield_data=yield_data(3:height(yield_data),:);


%Creating indices for different windows
t3 = datetime(1989,3,1);
t1 = datetime(1997,12,1);
t2 = datetime(2022,6,1);
t = t1:calmonths(1):t2;
t4 = datetime(2008,11,1);
t_pre_zlb = t1:calmonths(1):t4;
t5 = datetime(2015,11,1);
t_last = t5:calmonths(1):t2;
t_zlb=t4:calmonths(1):t5;

%Putting yield data in timetable format. Default window for yields: full
%sample
% Step 1: Convert 'yield_data.Var1' and 't' to datetime objects
data_dates = datetime(yield_data.Var1, 'InputFormat', 'yyyyMM');  % Correct conversion
target_dates = datetime(t, 'InputFormat', 'yyyyMM');              % Correct conversion

% Step 2: Strip time components to keep only year and month
data_dates = dateshift(data_dates, 'start', 'month');
target_dates = dateshift(target_dates, 'start', 'month');

% Step 3: Use ismember to find the matching dates
rowsToKeep = ismember(data_dates, target_dates);

yield_data = yield_data(rowsToKeep, :);
yield_data=table2timetable([table(datetime(t','Format','yyyyMM')), yield_data(:,2:size(yield_data,2))]);

%Now we will plot the yield data
figure;

% Extracting time vector
timeVector = yield_data.Var1;

% Assuming your timetable is named 'yourTimetable'
columnNames = yield_data.Properties.VariableNames;

% Plotting each yield column
hold on; % To overlay plots on the same figure
for i = 1:length(columnNames)
    hold on;  % Add 'hold on' to overlay plots on the same figure
    plot(timeVector, yield_data.(columnNames{i}), 'DisplayName', columnNames{i});
end
% plot(timeVector, yield_data.I02503MIndex, 'DisplayName', 'I02503MIndex');
% plot(timeVector, yield_data.I02506MIndex, 'DisplayName', 'I02506MIndex');
% Add more lines for other yields

hold off;

% Adding labels and legend
xlabel('Time');
ylabel('Yield');
title('Yield Curves Over Time');
legend('show');

% Extract maturity names and time indices
maturityNames = yield_data.Properties.VariableNames;
timeIndices = yield_data.Var1;

% Plot the yield residuals using surf
surf(table2array(yield_data));

% Set x-axis labels to be the maturity names
xticks(1:length(maturityNames));
xticklabels(maturityNames);

% Find the indices corresponding to the start of each year
yearStartIndices = find(month(timeIndices) == 1);

% Select every other year for the y-ticks
yearStartIndices = yearStartIndices(1:2:end);

% Get the corresponding years for labels (in datetime format)
ytickLabels = timeIndices(yearStartIndices);

% Convert datetime to cell array of strings for y-tick labels
ytickLabelsStr = cellstr(datestr(ytickLabels, 'yyyy'));

% Set the y-ticks to the start of each year
yticks(yearStartIndices);

% Set the y-tick labels to display only the years
yticklabels(ytickLabelsStr);

% Reverse the y-axis direction to go from past to recent
set(gca, 'YDir', 'reverse');

title('Yield Curve');
xlabel('Maturity');
ylabel('Time Index');
zlabel('Yields'); % Set the viewing angle for better visualization
%Now we will extract and plot Nelson-Siegel betas
[betas,ns_factor_loadings] = betas_function(yield_data);

%Plots of the betas
figure;

hold on
plot(betas.Var1,betas(:,1).betas1, 'DisplayName','Beta 1');
title('Beta 1')
plot(betas.Var1,betas(:,2).betas2, 'DisplayName','Beta 2');
title('Beta 2')
plot(betas.Var1,betas(:,3).betas3, 'DisplayName','Beta 3');
title('Beta 3')
hold off

% Adding labels and legend
xlabel('Time');
ylabel('Beta');
title('Betas Over Time');
legend('show');
%stationarity tests
[h,pValue] = adftest(betas(:,1));
[h,pValue] = adftest(betas(:,2));
[h,pValue] = adftest(betas(:,3));

[h,pValue] = kpsstest(betas(:,1));
[h,pValue] = kpsstest(betas(:,2));
[h,pValue] = kpsstest(betas(:,3));

[h,pValue] = pptest(betas(:,1));
[h,pValue] = pptest(betas(:,2));
[h,pValue] = pptest(betas(:,3));

%Now doing a table of the betas
table_betas=table_betas_function(betas)
writetable(table_betas,'betas_table.xlsx');


[estimated_yield_curve,yield_residuals]=yield_residuals_function(yield_data,ns_factor_loadings,betas);
surf(table2array(yield_residuals));
title('Yield Curve Residuals');

[estimated_yield_curve, yield_residuals] = yield_residuals_function(yield_data, ns_factor_loadings, betas);

% Extract maturity names and time indices
maturityNames = yield_data.Properties.VariableNames;
timeIndices = yield_residuals.Var1;

% Plot the yield residuals using surf
surf(table2array(yield_residuals));

% Set x-axis labels to be the maturity names
xticks(1:length(maturityNames));
xticklabels(maturityNames);

% Find the indices corresponding to the start of each year
yearStartIndices = find(month(timeIndices) == 1);

% Select every other year for the y-ticks
yearStartIndices = yearStartIndices(1:2:end);

% Get the corresponding years for labels (in datetime format)
ytickLabels = timeIndices(yearStartIndices);

% Convert datetime to cell array of strings for y-tick labels
ytickLabelsStr = cellstr(datestr(ytickLabels, 'yyyy'));

% Set the y-ticks to the start of each year
yticks(yearStartIndices);

% Set the y-tick labels to display only the years
yticklabels(ytickLabelsStr);

% Reverse the y-axis direction to go from past to recent
set(gca, 'YDir', 'reverse');

title('Yield Curve Residuals');
xlabel('Maturity');
ylabel('Time Index');
zlabel('Yield Residuals');

slope = table2array(yield_data(:,12))-table2array(yield_data(:,1));

nelson_siegel_slope = -0.7773*table2array(betas(:,2))+0.0551*table2array(betas(:,3));

corr(slope,nelson_siegel_slope);
%Now we will read the macro data
macro_data=readtable("stationary_data_for_macro_factors_2.csv");%Usada
%hasta 24-1-2025
% macro_data=readtable("stationary_data_for_macro_factors_3.csv");
macro_data=table2timetable([table(datetime(macro_data.Var1,'Format','yyyyMM')), macro_data(:,2:size(macro_data,2))]);

macro_data_stock=readtable("data_stock_and_watson_stationary_for_macro_factors_3.csv");
macro_data_stock=table2timetable([table(datetime(macro_data_stock.Var1,'Format','yyyyMM')), macro_data_stock(:,2:size(macro_data_stock,2))]);

macro_data_new_fast = readtable("new_fast_variables.csv");
macro_data_new_fast=table2timetable([table(datetime(macro_data_new_fast.Var1,'Format','yyyyMM')), macro_data_new_fast(:,2:size(macro_data_new_fast,2))]);

macro_data_final=synchronize(macro_data,macro_data_stock,macro_data_new_fast);
macro_data_final=rmmissing(macro_data_final);

macro_data_final.Federal_Surplus_or_Deficit____ = [NaN;100 * diff(macro_data_final.Federal_Surplus_or_Deficit____) ./ macro_data_final.Federal_Surplus_or_Deficit____(1:end-1)];

% TT.PctChange_VarName = 100 * diff(TT.VarName) ./ TT.VarName(1:end-1);
%Specify variables to exclude
exclude_vars = {'Persons_unemployed_15_to_26_weeks', 'Persons_unemployed_15__weeks', 'BAACorporateBondMinusGS10', 'x30_yearUSMortgageRateBondMinusGS10', 'USPMIManufacturing','cp_tb3m','PercentageChangeOfTotalNonrevolvingConsumerCredit_AtAnnualRate_','UnemploymentRate','vix'};  % Example variable names to exclude
% 
% Get variable names of the timetable
var_names = macro_data_final.Properties.VariableNames;

% Loop through each variable
for i = 1:numel(var_names)
    var_name = var_names{i};
    if ~ismember(var_name, exclude_vars)
        % Apply rolling sum of order 12
        macro_data_final.(var_name) = movsum(macro_data_final.(var_name), [11 0], 'omitnan');
    end
end

% Eliminate the first 11 rows
macro_data_final = macro_data_final(12:end, :);
%Readind the FFR data


% % Contraparte trimestral para obtener transformaciones intertrimestrales
% % Get variable names of the timetable
% var_names = macro_data_final.Properties.VariableNames;
% 
% % Loop through each variable
% for i = 1:numel(var_names)
%     var_name = var_names{i};
%     if ~ismember(var_name, exclude_vars)
% %         Apply rolling sum of order 3
%         macro_data_final.(var_name) = movsum(macro_data_final.(var_name), [2 0], 'omitnan');
%     end
% end
% 
% % Eliminate the first 11 rows
% macro_data_final = macro_data_final(3:end, :);

macro_data_final=synchronize(macro_data_final,betas);
macro_data_final=rmmissing(macro_data_final);
macro_data_final = macro_data_final(:, 1:end-3);

ffr=readtable("\ffr.csv");

% Extracting the date strings from the first column
dateStrings = ffr{:, 1};

% Converting date strings to datetime with 'MMM yyyy' format
dates = datetime(dateStrings, 'InputFormat', 'MMM yyyy');

% Converting datetime to the desired format 'yyyyMM'
datesFormatted = datetime(year(dates), month(dates), 1, 'Format', 'yyyyMM');

% Creating a timetable with the formatted dates and other columns
ffr = timetable(datesFormatted, ffr{:, 2:end});

%Reading ICM
icm=readtable("\Shadow_FFR_1221.xlsx");
icm=table2timetable([table(datetime(char(icm{:,1}),'Format','yyyyMM')), icm(:,2:size(icm,2))]);

% regression_variables = macro_data_final(:,[53 54 56])
% 
% icm_regression_data=synchronize(icm,regression_variables);
% icm_regression_data=rmmissing(icm_regression_data);
% 
% % Extract response and predictor variables correctly
% response_var = icm_regression_data{:, 1};
% predictors = icm_regression_data{:, [3, 4, 5]};
% 
% % Fit the linear model
% mdl = fitlm(predictors, response_var);
% 
% icm = mdl.Fitted;
% icm=table2timetable([table(datetime(icm_regression_data.Var1,'Format','yyyyMM')), array2table(icm)]);

%Now we will construct the Monetary Policy Indicator for the whole sample
values = ffr.Var1;

cutoffDate = datetime(icm.Var1(end));

belowThreshold = ffr(values < 0.25 & ffr.datesFormatted < cutoffDate, :);

ffr_icm = ffr;

for_dates = ffr_icm(values < 0.25 & ffr_icm.datesFormatted < cutoffDate, :);

valuesForDesiredDate = icm.Shadow_FFR(ismember(icm.Var1, for_dates.datesFormatted));

ffr_icm(values < 0.25 & ffr_icm.datesFormatted < cutoffDate, :) = array2table(valuesForDesiredDate);

%The whole data
%data=synchronize(macro_data,ffr_icm,betas);
%data=rmmissing(data);


slow_moving_variables=macro_data(:,2:47);
slow_moving_variables=[slow_moving_variables macro_data(:,59)];
slow_moving_variables=[slow_moving_variables macro_data(:,61:78)];

% slow_moving_variables_final=[macro_data_final(:,1:10) macro_data_final(:,12:18) macro_data_final(:,20:31) macro_data_final(:,33:47)];
% slow_moving_variables_final=[slow_moving_variables_final macro_data_final(:,59:82)];
% slow_moving_variables_final=[slow_moving_variables_final macro_data_final(:,85:93) macro_data_final(:,96:100) macro_data_final(:,103) macro_data_final(:,106)]; 

slow_moving_variables_final=[macro_data_final(:,1:10) macro_data_final(:,12:18) macro_data_final(:,20:31) macro_data_final(:,33:47)];
slow_moving_variables_final=[slow_moving_variables_final macro_data_final(:,59:82)];
slow_moving_variables_final=[slow_moving_variables_final macro_data_final(:,85:93) macro_data_final(:,103)]; 

% slow_moving_variables_final=[macro_data_final(:,1:10) macro_data_final(:,12:18) macro_data_final(:,20:31) macro_data_final(:,33:40)];
% slow_moving_variables_final=[slow_moving_variables_final macro_data_final(:,59:82)];
% slow_moving_variables_final=[slow_moving_variables_final macro_data_final(:,85:93) macro_data_final(:,103)]; 

%Now we do preliminary factor analysis for the slow moving variables
[coeff,score,latent,tsquared,explained,mu] = pca(zscore(table2array(slow_moving_variables_final)));
slow_macro_factors=zscore(table2array(slow_moving_variables_final))*coeff;
number_of_slow_factors=2;
slow_macro_factors=slow_macro_factors(:,1:number_of_slow_factors);
slow_macro_factors=array2timetable(array2table(slow_macro_factors),'RowTimes', datetime(macro_data_final.Var1,'Format','yyyyMM'));
std(table2array(table2array(slow_macro_factors)))


bar(explained)
ylabel('Marginal R2')
xlabel('Number of factors')
legend('Marginal R2 for slow-moving variables', 'Location', 'northeast')

plot(slow_macro_factors.Time,table2array(slow_macro_factors.Var1), 'DisplayName', 'First Slow Macro Factor');
hold on 
plot(slow_macro_factors.Time,table2array(slow_macro_factors.Var2), 'DisplayName', 'Second Slow Macro Factor');

legend('show');


fast_moving_variables=macro_data(:,48:58);

% fast_moving_variables_final=macro_data_final(:,48:54);
% fast_moving_variables_final=[fast_moving_variables_final macro_data_final(:,56:58) macro_data_final(:,83:84) macro_data_final(:,95) macro_data_final(:,101:102) macro_data_final(:,104:105) macro_data_final(:,107) macro_data_final(:,108:end)]; 

fast_moving_variables_final=macro_data_final(:,48:54);
fast_moving_variables_final=[fast_moving_variables_final macro_data_final(:,56:57) macro_data_final(:,83:84) macro_data_final(:,95) macro_data_final(:,101:102) macro_data_final(:,104:105) macro_data_final(:,107) macro_data_final(:,108:end) macro_data_final(:,96:100) macro_data_final(:,106)]; %macro_data_final(:,56:58)

% fast_moving_variables_final=macro_data_final(:,48:54);
% fast_moving_variables_final=[fast_moving_variables_final macro_data_final(:,56:58) macro_data_final(:,83:84) macro_data_final(:,95) macro_data_final(:,101:102) macro_data_final(:,104:105) macro_data_final(:,107) macro_data_final(:,108:end) macro_data_final(:,96:100) macro_data_final(:,106) macro_data_final(:,41:47)]; 

%Now we do preliminary factor analysis for the fast moving variables
[coeff,score,latent,tsquared,explained,mu] = pca(zscore(table2array(fast_moving_variables_final)));
fast_macro_factors=zscore(table2array(fast_moving_variables_final))*coeff;
number_of_fast_macro_factors=2;
fast_macro_factors=fast_macro_factors(:,1:number_of_fast_macro_factors);
fast_macro_factors=array2timetable(array2table(fast_macro_factors),'RowTimes', datetime(macro_data_final.Var1,'Format','yyyyMM'));
std(table2array(table2array(fast_macro_factors)))

bar(explained)
ylabel('Marginal R2')
xlabel('Number of factors')
legend('Marginal R2 for fast-moving variables', 'Location', 'northeast')

plot(fast_macro_factors.Time,table2array(fast_macro_factors.Var1), 'DisplayName', 'First Fast Macro Factor');
hold on 
plot(fast_macro_factors.Time,table2array(fast_macro_factors.Var2), 'DisplayName', 'Second Fast Macro Factor');

legend('show');


%let's reorder slow moving variables to identify the shocks correctly with
%the name factor normalization
%At the end we took off a sectorial IP and a sectorial CPI to avoid
%multicollinearity problems
slow_moving_variables_final=[slow_moving_variables_final(:,1) slow_moving_variables_final(:,47) slow_moving_variables_final(:,2) slow_moving_variables_final(:,4:46) slow_moving_variables_final(:,48:50) slow_moving_variables_final(:,52:end)];
%reordering fast moving variables
fast_moving_variables_final=[fast_moving_variables_final(:,1) fast_moving_variables_final(:,4) fast_moving_variables_final(:,2:3) fast_moving_variables_final(:,5:end)];
matrix_inverted_slow=inv(table2array(table2array(slow_macro_factors))'*table2array(table2array(slow_macro_factors)));

%Now we set up the dataset for model estimation
model_estimation_variables=synchronize(slow_moving_variables_final,ffr_icm,fast_moving_variables_final,betas);
% Extract the variables (columns) as an array
dataArray = betas{:, :};

% Calculate the interannual difference by shifting by 12 months
interannualDiffArray = dataArray(13:end, :) - dataArray(1:end-12, :);

% Create a new timetable with the adjusted time vector
interannualDiffbetas = array2timetable(interannualDiffArray, 'RowTimes', betas.Var1(13:end), 'VariableNames', betas.Properties.VariableNames);

dataArray = ffr_icm{:, 1};

% Calculate the interannual difference by shifting by 12 months
interannualDiffArray = dataArray(13:end) - dataArray(1:end-12);

% Create a new timetable with the adjusted time vector
interannualDiffffr_icm = array2timetable(interannualDiffArray, 'RowTimes', ffr_icm.datesFormatted(13:end), 'VariableNames', ffr_icm.Properties.VariableNames);

% Extract the data
dataArray = betas{:, :};

% Calculate the differences
diffArray = diff(dataArray);

% Adjust the time vector to match the result of the differences
newTime = betas.Var1(2:end);

% Create a new timetable with the differences
diffbetas = array2timetable(diffArray, 'RowTimes', newTime, 'VariableNames', betas.Properties.VariableNames);

% Extract the data
dataArray = ffr_icm{:, :};

% Calculate the differences
diffArray = diff(dataArray);

% Adjust the time vector to match the result of the differences
newTime = ffr_icm.datesFormatted(2:end);

% Create a new timetable with the differences
diffffr_icm = array2timetable(diffArray, 'RowTimes', newTime, 'VariableNames', ffr_icm.Properties.VariableNames);

%Let's detrend beta 1

% %Extract beta 1 as a double
% beta_1 = betas.betas1;
% 
% %Generate time as a double starting in t=1
% time = (1:size(betas,1))'; 
% 
% % Fit a linear trend
% p = polyfit(time, beta_1, 1); % Linear fit: p(1) = slope, p(2) = intercept
% 
% % Compute the trend
% trend = polyval(p, time);
% 
% % Detrend the series
% detrended_beta_1 = beta_1 - trend;
% 
% % Add the detrended series back to the timetable
% betas.betas1 = detrended_beta_1;


% model_estimation_variables=synchronize(slow_moving_variables_final,ffr_icm,fast_moving_variables_final,betas);
% model_estimation_variables=rmmissing(model_estimation_variables);
% %Putting the variables in the right order for structural model estimation
% model_estimation_variables = model_estimation_variables(:, [3,2, 1,4:end]);

slow_moving_variables_final_cleaned = varfun(@(data) ...
    replaceOutliersWithMovingAverage(data), ...
    slow_moving_variables_final, 'InputVariables', slow_moving_variables_final.Properties.VariableNames);

fast_moving_variables_final_cleaned = varfun(@(data) ...
    replaceOutliersWithMovingAverage(data), ...
    fast_moving_variables_final, 'InputVariables', fast_moving_variables_final.Properties.VariableNames);

model_estimation_variables=synchronize(slow_moving_variables_final_cleaned,ffr_icm,fast_moving_variables_final_cleaned,betas);
model_estimation_variables=rmmissing(model_estimation_variables);
%Putting the variables in the right order for structural model estimation
model_estimation_variables = model_estimation_variables(:, [3,2, 1,4:end]);
%Example parameters
est_par.smpl_par.calvec = 1:size(model_estimation_variables,1);
est_par.smpl_par.nper = 12;  % Monthly data
est_par.smpl_par.nfirst = [1,1]; % Index of the first observation
est_par.smpl_par.nlast = [length(est_par.smpl_par.calvec),length(est_par.smpl_par.calvec)];  % Last observation for estimation

est_par.fac_par.w = [model_estimation_variables.Var1,model_estimation_variables.betas1,model_estimation_variables.betas2,model_estimation_variables.betas3];  % Matrix of observed factors (if applicable)
est_par.fac_par.lambda_constraints= 1;  % Constraint matrix
est_par.fac_par.nt_min = 50;  % Minimum number of observations for factor estimation
est_par.fac_par.tol = 1e-6;  % Precision of the estimate
est_par.fac_par.nfac.unobserved = 4;  % Number of unobserved factors
est_par.fac_par.nfac.observed = 4;  % Number of observed factors
est_par.fac_par.nfac.total = est_par.fac_par.nfac.unobserved + est_par.fac_par.nfac.observed;
est_par.n_uarlag = 1; %number of AR lags for uniquesess
est_par.lambda.nt_min = 50;
est_par.fac_par.lambda_constraints_full = 1;

n_naming_slow_variables = 2;

n_non_naming_slow_variables = size(slow_moving_variables_final,2)-n_naming_slow_variables;%80;
% Initialize the matrix
matrix_restrictions_slow = zeros(n_non_naming_slow_variables, est_par.fac_par.nfac.total+2);

% Fill the upper left 6x6 block with the desired pattern
matrix_restrictions_slow (1:6,:) = [3, 1, 0, 0, 0, 0, 0, 0, 0, 0 ;
    3, 0, 1, 0, 0, 0, 0, 0, 0, 0 ;
    3, 0, 0, 1, 0, 0, 0, 0, 0, 0 ;
    3, 0, 0, 0, 1, 0, 0, 0, 0, 0 ;
    3, 0, 0, 0, 0, 0, 0, 1, 0, 0 ;
    3, 0, 0, 0, 0, 0, 0, 0, 1, 0 ];

% Replace the number 3 with different values (from 4 to 65)
for j = 4:(n_non_naming_slow_variables+2)
    % Copy the pattern
    pattern = matrix_restrictions_slow (1:6, :);
    
    % Replace 3 with the current value of j
    pattern(pattern == 3) = j;
    
    % Add the modified pattern to the matrix
    matrix_restrictions_slow (((j-3)*6+1):(j-2)*6,:) = pattern;
end

% % Fill the upper left 6x6 block with the desired pattern
% matrix_restrictions_slow (1:6,:) = [5, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  ;
%     5, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 ;
%     5, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 ;
%     5, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 ;
%     5, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 ;
%     5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 ];
% 
% % Replace the number 3 with different values (from 4 to 65)
% for j = 6:65
%     % Copy the pattern
%     pattern = matrix_restrictions_slow (1:6, :);
%     
%     % Replace 3 with the current value of j
%     pattern(pattern == 5) = j;
%     
%     % Add the modified pattern to the matrix
%     matrix_restrictions_slow (((j-5)*6+1):(j-4)*6,:) = pattern;
% end




% Replace the number 3 with different values (from 4 to 65)
% for j = 69:76
%     % Copy the pattern
%     pattern = matrix_restrictions_fast (1:3, :);
%     
%     % Replace 3 with the current value of j
%     pattern(pattern == 68) = j;
%     
%     % Add the modified pattern to the matrix
%     matrix_restrictions_fast (((j-68)*3+1):(j-67)*3,:) = pattern;
% end

est_par.fac_par.lambda_constraints_est = 1;
est_par.fac_par.lambda_constraints_est = [
    1, 1, 0, 0, 0, 0, 0, 0, 0, 0 ;
    1, 0, 1, 0, 0, 0, 0, 0, 0, 0 ;
    1, 0, 0, 1, 0, 0, 0, 0, 0, 0 ;
    1, 0, 0, 0, 1, 0, 0, 0, 0, 0 ;
    1, 0, 0, 0, 0, 1, 0, 0, 0, 1;
    1, 0, 0, 0, 0, 0, 1, 0, 0, 0;
    1, 0, 0, 0, 0, 0, 0, 1, 0, 0 ;
    1, 0, 0, 0, 0, 0, 0, 0, 1, 0 ;
    2, 1, 0, 0, 0, 0, 0, 0, 0, 0 ;
    2, 0, 1, 0, 0, 0, 0, 0, 0, 0 ;
    2, 0, 0, 1, 0, 0, 0, 0, 0, 0 ;
    2, 0, 0, 0, 1, 0, 0, 0, 0, 0 ;
    2, 0, 0, 0, 0, 1, 0, 0, 0, 0;
    2, 0, 0, 0, 0, 0, 1, 0, 0, 1;
    2, 0, 0, 0, 0, 0, 0, 1, 0, 0 ;
    2, 0, 0, 0, 0, 0, 0, 0, 1, 0 ;
    matrix_restrictions_slow ;
    83, 0, 1, 0, 0, 0, 0, 0, 0, 0 ;
    83, 0, 0, 1, 0, 0, 0, 0, 0, 0 ;
    83, 0, 0, 0, 1, 0, 0, 0, 0, 0 ;
    83, 0, 0, 0, 0, 1, 0, 0, 0, 0;
    83, 0, 0, 0, 0, 0, 1, 0, 0, 0;
    83, 0, 0, 0, 0, 0, 0, 1, 0, 1;
    83, 0, 0, 0, 0, 0, 0, 0, 1, 0 ;
    84, 0, 1, 0, 0, 0, 0, 0, 0, 0 ;
    84, 0, 0, 1, 0, 0, 0, 0, 0, 0 ;
    84, 0, 0, 0, 1, 0, 0, 0, 0, 0 ;
    84, 0, 0, 0, 0, 1, 0, 0, 0, 0;
    84, 0, 0, 0, 0, 0, 1, 0, 0, 0;
    84, 0, 0, 0, 0, 0, 0, 1, 0, 0 ;
    84, 0, 0, 0, 0, 0, 0, 0, 1, 1 ;
    ]%matrix_restrictions_fast

% est_par.fac_par.lambda_constraints_est = [
%     1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ;
%     1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0  ;
%     1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0  ;
%     1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0  ;
%     1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 , 1;
%     1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 ;
%     1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0  ;
%     1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0  ;
%     1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0  ;
%     1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0  ;
%     2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  ;
%     2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0  ;
%     2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0  ;
%     2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0  ;
%     2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 ;
%     2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 , 1;
%     2, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0  ;
%     2, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0  ;
%     2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0  ;
%     2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0  ;
%     3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  ;
%     3, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0  ;
%     3, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0  ;
%     3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0  ;
%     3, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 ;
%     3, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 , 0;
%     3, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1  ;
%     3, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0  ;
%     3, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0  ;
%     3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0  ;
%     4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  ;
%     4, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0  ;
%     4, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0  ;
%     4, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0  ;
%     4, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 ;
%     4, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 , 0;
%     4, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0  ;
%     4, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1  ;
%     4, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0  ;
%     4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0  ;
%     matrix_restrictions_slow ;
%     66, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
%     66, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
%     66, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
%     66, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 ;
%     66, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
%     66, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
%     66, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
%     66, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 ;
%     66, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1;
%     66, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 ;
%     67, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
%     67, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
%     67, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
%     67, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 ;
%     67, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;
%     67, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 ]
est_par.fac_par.lambda_constraints_full = [
    1, 1, 0, 0, 0, 0, 0, 0, 0, 0 ;
    1, 0, 1, 0, 0, 0, 0, 0, 0, 0 ;
    1, 0, 0, 1, 0, 0, 0, 0, 0, 0 ;
    1, 0, 0, 0, 1, 0, 0, 0, 0, 0 ;
    1, 0, 0, 0, 0, 1, 0, 0, 0, 1;
    1, 0, 0, 0, 0, 0, 1, 0, 0, 0;
    1, 0, 0, 0, 0, 0, 0, 1, 0, 0 ;
    1, 0, 0, 0, 0, 0, 0, 0, 1, 0 ;
    2, 1, 0, 0, 0, 0, 0, 0, 0, 0 ;
    2, 0, 1, 0, 0, 0, 0, 0, 0, 0 ;
    2, 0, 0, 1, 0, 0, 0, 0, 0, 0 ;
    2, 0, 0, 0, 1, 0, 0, 0, 0, 0 ;
    2, 0, 0, 0, 0, 1, 0, 0, 0, 0;
    2, 0, 0, 0, 0, 0, 1, 0, 0, 1;
    2, 0, 0, 0, 0, 0, 0, 1, 0, 0 ;
    2, 0, 0, 0, 0, 0, 0, 0, 1, 0 ;
    matrix_restrictions_slow;
    83, 0, 1, 0, 0, 0, 0, 0, 0, 0 ;
    83, 0, 0, 1, 0, 0, 0, 0, 0, 0 ;
    83, 0, 0, 0, 1, 0, 0, 0, 0, 0 ;
    83, 0, 0, 0, 0, 1, 0, 0, 0, 0;
    83, 0, 0, 0, 0, 0, 1, 0, 0, 0;
    83, 0, 0, 0, 0, 0, 0, 1, 0, 1;
    83, 0, 0, 0, 0, 0, 0, 0, 1, 0 ;
    84, 0, 1, 0, 0, 0, 0, 0, 0, 0 ;
    84, 0, 0, 1, 0, 0, 0, 0, 0, 0 ;
    84, 0, 0, 0, 1, 0, 0, 0, 0, 0 ;
    84, 0, 0, 0, 0, 1, 0, 0, 0, 0;
    84, 0, 0, 0, 0, 0, 1, 0, 0, 0;
    84, 0, 0, 0, 0, 0, 0, 1, 0, 0 ;
    84, 0, 0, 0, 0, 0, 0, 0, 1, 1;
    ]%matrix_restrictions_fast

% est_par.fac_par.lambda_constraints_full = [
%     1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ;
%     1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0  ;
%     1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0  ;
%     1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0  ;
%     1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 , 1;
%     1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 ;
%     1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0  ;
%     1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0  ;
%     1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0  ;
%     1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0  ;
%     2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  ;
%     2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0  ;
%     2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0  ;
%     2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0  ;
%     2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 ;
%     2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 , 1;
%     2, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0  ;
%     2, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0  ;
%     2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0  ;
%     2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0  ;
%     3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  ;
%     3, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0  ;
%     3, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0  ;
%     3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0  ;
%     3, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 ;
%     3, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 , 0;
%     3, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1  ;
%     3, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0  ;
%     3, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0  ;
%     3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0  ;
%     4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  ;
%     4, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0  ;
%     4, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0  ;
%     4, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0  ;
%     4, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 ;
%     4, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 , 0;
%     4, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0  ;
%     4, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1  ;
%     4, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0  ;
%     4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0  ;
%     matrix_restrictions_slow ;
%     66, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
%     66, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
%     66, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
%     66, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 ;
%     66, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
%     66, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
%     66, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
%     66, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 ;
%     66, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1;
%     66, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 ;
%     67, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
%     67, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
%     67, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
%     67, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 ;
%     67, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;
%     67, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 ]
% columnsToKeep = setdiff(model_estimation_variables.Properties.VariableNames, {'Var1'});
% subsetTimetable = model_estimation_variables(:, columnsToKeep);
% est_data = timetable2table(subsetTimetable);
% est_data = table2array(est_data(:, 2:end));

% Find the index of 'Var1' in the variable names
idxVar1 = strcmp(model_estimation_variables.Properties.VariableNames, 'Var1');

% Get the names of the variables excluding 'Var1'
columnsToKeep = model_estimation_variables.Properties.VariableNames(~idxVar1);

% Extract the variables in the same order as they appear in the original table
variablesToKeep = model_estimation_variables(:, columnsToKeep);

est_data = timetable2table(variablesToKeep);
est_data = table2array(est_data(:, 2:end));

% Names of variables to exclude
excludeVars = {'Var1', 'betas1', 'betas2', 'betas3'};

% Get the indices of variables to exclude
idxExclude = ismember(model_estimation_variables.Properties.VariableNames, excludeVars);

% Get the names of variables to keep
columnsToKeep = model_estimation_variables.Properties.VariableNames(~idxExclude);

% Extract the variables to keep
variablesToKeep = model_estimation_variables(:, ~idxExclude);

est_data = timetable2table(variablesToKeep);
est_data = table2array(est_data(:, 2:end));


lsout = factor_estimation_ls(est_data, est_par);


bn_icp = bai_ng(lsout);

data=est_data;

%est_par.fac_par.lag_order = 5;
var_par.nlag   = 2;
var_par.icomp  = 1;
var_par.iconst = 1;
est_par.var_par.nlag   = 2;
est_par.var_par.icomp  = 1;
est_par.var_par.iconst = 1;
out = factor_estimation_ls_full(data, est_par, var_par);

%Diagnosis 
n_of_factors = size(out.fac,2);
resid = out.varout.resid;
for i = 1:n_of_factors
table_residuals_autocorrelation_test(i,:) = lbqtest(resid((var_par.nlag+1):end,i));
end
table_residuals_autocorrelation_test;
multivariate_test_lags = 4;
[h,pValue,stat,cValue]=mlbqtest(resid((var_par.nlag+1):end,:),multivariate_test_lags);

AICs=lag_order_selection(out.fac,6);
[min_num_AIC,min_idx_AIC] = min(AICs);

BICs=lag_order_selection_BIC(out.fac,6);
[min_num_BIC,min_idx_BIC] = min(BICs);

%Up until here well
factors = out.fac;

% Assuming your 365x3 double array is named 'yourData'
data = factors;

% Assuming your date column is named 'DateColumn'
dateColumn = model_estimation_variables.Var1_1; % Replace 'yourTable' with your actual table variable

% Convert the date column to datetime
rowTimes = datetime(dateColumn, 'InputFormat', 'yyyyMM');

% Create default variable names ('Var1', 'Var2', 'Var3')
variableNames = {'Factor1', 'Factor2', 'Factor3', 'Factor4'};

% Create timetable
factors = timetable(rowTimes, data(:, 1), data(:, 2), data(:, 3), data(:, 4), 'VariableNames', variableNames);

% Plotting the three time series on the same graph
figure;
plot(factors.rowTimes,factors.Factor1, '-o', 'DisplayName', 'Factor1');
hold on;
plot(factors.rowTimes,factors.Factor2, '-s', 'DisplayName', 'Factor2');
hold on;
plot(factors.rowTimes,factors.Factor3, '-d', 'DisplayName', 'Factor3');
hold on;
plot(factors.rowTimes,factors.Factor4, '-d', 'DisplayName', 'Factor4');

% Add labels and title
xlabel('Time');
ylabel('Value');
title('Factors Time Series Plot');
legend('show');


[h,pValue] = adftest(factors.Factor1);

[h,pValue] = kpsstest(factors.Factor1);

[h,pValue] = pptest(factors.Factor1);

lag_order = 2;
n_factors = 8;

fac_est_out = out;
fac_est_out.fac = fac_est_out.fac(:, [5:6,1,7:8,2:4]);%fac_est_out.fac(:, [5:6,1,7:8,2:4]);
fac_est_out.lam_mat = fac_est_out.lam_mat(:, [5:6,1,7:8,2:4]);%fac_est_out.lam_mat(:, [5:6,1,7:8,2:4]);
fac_est_out.varout.betahat = fac_est_out.varout.betahat(:, [5:6,1,7:8,2:4]);%fac_est_out.varout.betahat(:, [5:6,1,7:8,2:4]);
fac_est_out.varout.seps = fac_est_out.varout.seps(:, [5:6,1,7:8,2:4]);%fac_est_out.varout.seps(:, [5:6,1,7:8,2:4]);
fac_est_out.varout.seps = fac_est_out.varout.seps([5:6,1,7:8,2:4], :);
fac_est_out.varout.resid = fac_est_out.varout.resid(:, [5:6,1,7:8,2:4]);
% fac_est_out.varout.coef.Q = fac_est_out.varout.coef.Q([5:8,1:4],: );
fac_est_out.varout.coef.M(1:8,:) = fac_est_out.varout.coef.M(1:8,reshape([[5:6,1,7:8,2:4]+[0 + (0:(lag_order-1)) *8]']', 1, []));%fac_est_out.varout.coef.M(1:8,[5:6,1,7:8,2:4,13,14,9,15,16,10,11,12,21,22,17,23,24,18,19,20]);
fac_est_out.varout.coef.M= fac_est_out.varout.coef.M([5:6,1,7:8,2:4,9:(lag_order*n_factors)],:);%fac_est_out.varout.coef.M([5:6,1,7:8,2:4,9:24],:);
fac_est_out.varout.coef.G = [chol(fac_est_out.varout.seps).';fac_est_out.varout.coef.G(9:(lag_order*n_factors),:)];%[chol(fac_est_out.varout.seps).';fac_est_out.varout.coef.G(9:24,:)];

fac_est_par = est_par;
fac_est_par.fac_par.lambda_constraints_est = fac_est_par.fac_par.lambda_constraints_est(:, [5:6,1,7:8,2:4]+1);%fac_est_par.fac_par.lambda_constraints_est(:, [6:7,2,8:9,3:5]);
fac_est_par.fac_par.lambda_constraints_full = fac_est_par.fac_par.lambda_constraints_full(:, [5:6,1,7:8,2:4]+1);%fac_est_par.fac_par.lambda_constraints_full(:, [6:7,2,8:9,3:5]);

% fac_est_out = out;
% fac_est_out.fac = fac_est_out.fac(:, [5:8,1,9:10,2:4]);
% fac_est_out.lam_mat = fac_est_out.lam_mat(:, [5:8,1,9:10,2:4]);
% fac_est_out.varout.betahat = fac_est_out.varout.betahat(:, [5:8,1,9:10,2:4]);
% fac_est_out.varout.seps = fac_est_out.varout.seps(:, [5:8,1,9:10,2:4]);
% fac_est_out.varout.seps = fac_est_out.varout.seps([5:8,1,9:10,2:4], :);
% fac_est_out.varout.resid = fac_est_out.varout.resid(:, [5:8,1,9:10,2:4]);
% % fac_est_out.varout.coef.Q = fac_est_out.varout.coef.Q([5:8,1:4],: );
% fac_est_out.varout.coef.M(1:10,:) = fac_est_out.varout.coef.M(1:10,[5:8,1,9:10,2:4,15:18,11,19:20,12:14,25:28,21,29:30,22:24]);%fac_est_out.varout.coef.M(1:8,[5:6,1,7:8,2:4,13,14,9,15,16,10,11,12,21,22,17,23,24,18,19,20,29,30,25,31,32,26,27,28,37,38,33,39,40,34,35,36,45,46,41,47,48,42,43,44,53,54,49,55,56,50,51,52,61,62,57,63,64,58,59,60])
% fac_est_out.varout.coef.M= fac_est_out.varout.coef.M([5:8,1,9:10,2:4,11:30],:);
% fac_est_out.varout.coef.G = [chol(fac_est_out.varout.seps).';fac_est_out.varout.coef.G(11:30,:)];
% 
% fac_est_par = est_par;
% fac_est_par.fac_par.lambda_constraints_est = fac_est_par.fac_par.lambda_constraints_est(:, [6:9,2,10:11,3:5])
% fac_est_par.fac_par.lambda_constraints_full = fac_est_par.fac_par.lambda_constraints_full(:, [6:9,2,10:11,3:5])

plot(fac_est_out.fac(:,5))

decomp_par.hor = 60;
decomp_par.varcum = 0;
decomp_par.cancor = 0;

%Number of observable variables panel
n_variables = size(est_data,2);%91 %100;

bptcodevec = ones(1, n_variables);%bptcodevec = ones(1, 97)

irf_vdecomp_out = dynamic_factor_irf_vdecomp(fac_est_out, est_par, decomp_par, bptcodevec);

%irf_vdecomp_out_factors = dynamic_factor_irf_vdecomp_factors(fac_est_out, est_par, decomp_par, bptcodevec);

% Extracting the data and reshaping it
data_irf = squeeze(irf_vdecomp_out.imp_y_fac_mat_scl(1,5, :));

% Plotting the reshaped data
plot(data_irf);

% -- Compute Standard Errors for IRFs and VDs using parametric bootstrap simulations
n_rep = 500;   % Number of bootstap simulations for computing SEs
datain.bpdata = est_data;
%n_variables = 97;
datain.bpinclcode = ones(1, n_variables);
datain.bptcodevec = ones(1, n_variables);
se_irf_vdecomp_out = se_dynamic_factor_irf_vdecomp(datain,out,est_par,decomp_par,n_rep,var_par);%se_dynamic_factor_irf_vdecomp(datain,out,est_par,decomp_par,n_rep,var_par);
%IRFs with restrictions
fac_est_out_irf_restrictions = fac_est_out;
%fac_est_out_irf_restrictions.varout.betahat([4,12,20],[6,7,8]) = 0;
fac_est_out_irf_restrictions.varout.coef.M([6,7,8],[3+(0:(lag_order-1)) *8]) = 0;%fac_est_out_irf_restrictions.varout.coef.M([6,7,8],[3,11,19])
%fac_est_out_irf_restrictions.varout.coef.M = fac_est_out_irf_restrictions.varout.coef.M*0.95
%fac_est_out_irf_restrictions.varout.coef.G(1,:) = 0;
%fac_est_out_irf_restrictions.varout.seps(1,:) = 0;

out_irf_restrictions = out;
%fac_est_out_irf_restrictions.varout.betahat([4,12,20],[6,7,8]) = 0;%[1,9,17]
out_irf_restrictions.varout.coef.M([2,3,4],[1+(0:(lag_order-1)) *8]) = 0;%out_irf_restrictions.varout.coef.M([2,3,4],[1,9,17]) 
%out_irf_restrictions.varout.coef.M = out_irf_restrictions.varout.coef.M*0.95
%fac_est_out_irf_restrictions.varout.coef.G(1,:) = 0;
%fac_est_out_irf_restrictions.varout.seps(1,:) = 0;

irf_vdecomp_out_restrictions = dynamic_factor_irf_vdecomp(fac_est_out_irf_restrictions, est_par, decomp_par, bptcodevec);

se_irf_vdecomp_out_restrictions = se_dynamic_factor_irf_vdecomp(datain,out_irf_restrictions,est_par,decomp_par,n_rep,var_par);%se_dynamic_factor_irf_vdecomp(datain,out,est_par,decomp_par,n_rep,var_par);

data_irf = squeeze(irf_vdecomp_out.imp_y_fac_mat_scl(1,5, :));

data_irf_restrictions = squeeze(irf_vdecomp_out_restrictions.imp_y_fac_mat_scl(1,5, :));%End of IRFs with restrictions

% Plotting the reshaped data
plot(data_irf_restrictions);



pointEstimate = plot(0.25*squeeze(irf_vdecomp_out.imp_y_fac_mat_scl(2,3, :)));
hold on
upper1 =        plot(0.25*squeeze(irf_vdecomp_out.imp_y_fac_mat_scl(2,3, :))+0.25*squeeze(se_irf_vdecomp_out.se_imp_y_fac_mat_scl(2,3, :))*1.5, '-y');
hold on
lower1 =        plot(0.25*squeeze(irf_vdecomp_out.imp_y_fac_mat_scl(2,3, :))-0.25*squeeze(se_irf_vdecomp_out.se_imp_y_fac_mat_scl(2,3, :))*1.5, '-y');
hold on
upper2 =        plot(0.25*squeeze(irf_vdecomp_out.imp_y_fac_mat_scl(2,3, :))+0.25*squeeze(se_irf_vdecomp_out.se_imp_y_fac_mat_scl(2,3, :))*2, '-r');
hold on
lower2 =        plot(0.25*squeeze(irf_vdecomp_out.imp_y_fac_mat_scl(2,3, :))-0.25*squeeze(se_irf_vdecomp_out.se_imp_y_fac_mat_scl(2,3, :))*2, '-r');
%hold on
%plot(0.25*squeeze(irf_vdecomp_out_restrictions.imp_y_fac_mat_scl(92,3, :)));
title('Impulse Response Function: Inflation')
xlabel('Time period')
ylabel('Effect')
legend([pointEstimate, upper1, upper2], ...
       'Point Estimate', 'Confidence Bands: 1.5 standard errors', 'Confidence Bands: 2 standard errors', ...
       'Location', 'Best');
yline(0,'--black','HandleVisibility','off')

pointEstimate = plot(0.25*squeeze(irf_vdecomp_out_restrictions.imp_y_fac_mat_scl(15,3, :)));
hold on
upper1 =        plot(0.25*squeeze(irf_vdecomp_out_restrictions.imp_y_fac_mat_scl(15,3, :))+0.25*squeeze(se_irf_vdecomp_out_restrictions.se_imp_y_fac_mat_scl(15,3, :))*1.5, '-y');
hold on
lower1 =        plot(0.25*squeeze(irf_vdecomp_out_restrictions.imp_y_fac_mat_scl(15,3, :))-0.25*squeeze(se_irf_vdecomp_out_restrictions.se_imp_y_fac_mat_scl(15,3, :))*1.5, '-y');
hold on
upper2 =        plot(0.25*squeeze(irf_vdecomp_out_restrictions.imp_y_fac_mat_scl(15,3, :))+0.25*squeeze(se_irf_vdecomp_out_restrictions.se_imp_y_fac_mat_scl(15,3, :))*2, '-r');
hold on
lower2 =        plot(0.25*squeeze(irf_vdecomp_out_restrictions.imp_y_fac_mat_scl(15,3, :))-0.25*squeeze(se_irf_vdecomp_out_restrictions.se_imp_y_fac_mat_scl(15,3, :))*2, '-r');
%hold on
%plot(0.25*squeeze(irf_vdecomp_out_restrictions.imp_y_fac_mat_scl(92,3, :)));
title('Impulse Response Function: Capacity utilization, total')
xlabel('Time period')
ylabel('Effect')
legend([pointEstimate, upper1, upper2], ...
       'Point Estimate', 'Confidence Bands: 1.5 standard errors', 'Confidence Bands: 2 standard errors', ...
       'Location', 'Best');
yline(0,'--black','HandleVisibility','off')


irf_factors = dynamic_factor_irf_vdecomp_factors(fac_est_out, est_par, decomp_par, bptcodevec);

se_irf_factors = se_dynamic_factor_irf_vdecomp_factors(datain,out,est_par,decomp_par,n_rep,var_par);%se_dynamic_factor_irf_vdecomp_factors(datain,out,est_par,decomp_par,n_rep,var_par);

data_factors_irf = squeeze(irf_factors.imp_y_fac_mat_scl(5,5, :));

% Plotting the reshaped data
plot(data_factors_irf);

pointEstimate = plot(-0.25*squeeze(irf_factors.imp_y_fac_mat_scl(7,3, :)));
hold on
upper1 =        plot(-0.25*squeeze(irf_factors.imp_y_fac_mat_scl(7,3, :))+0.25*squeeze(se_irf_factors.se_imp_y_fac_mat_scl(7,3, :))*1.5,'-y');
hold on
lower1 =        plot(-0.25*squeeze(irf_factors.imp_y_fac_mat_scl(7,3, :))-0.25*squeeze(se_irf_factors.se_imp_y_fac_mat_scl(7,3, :))*1.5,'-y');
hold on
upper2 =        plot(-0.25*squeeze(irf_factors.imp_y_fac_mat_scl(7,3, :))+0.25*squeeze(se_irf_factors.se_imp_y_fac_mat_scl(7,3, :))*2,'-r');
hold on
lower2 =        plot(-0.25*squeeze(irf_factors.imp_y_fac_mat_scl(7,3, :))-0.25*squeeze(se_irf_factors.se_imp_y_fac_mat_scl(7,3, :))*2,'-r');
title('Impulse Response Function: Yield Curve Slope (Negative of Nelson-Siegel Beta 2)')
xlabel('Time period')
ylabel('Effect')
legend([pointEstimate, upper1, upper2], ...
       'Point Estimate', 'Confidence Bands: 1.5 standard errors', 'Confidence Bands: 2 standard errors', ...
       'Location', 'Best');
yline(0,'--black','HandleVisibility','off')

plot(0.25*squeeze(irf_factors.imp_y_fac_mat_scl(6,3, :)+0.1367*irf_factors.imp_y_fac_mat_scl(7,3, :)+0.1361*irf_factors.imp_y_fac_mat_scl(8,3, :)-(irf_factors.imp_y_fac_mat_scl(6,3, :)+0.9140*irf_factors.imp_y_fac_mat_scl(7,3, :)+0.0810*irf_factors.imp_y_fac_mat_scl(8,3, :))));
hold on
plot(0.25*squeeze(irf_factors.imp_y_fac_mat_scl(6,3, :)+0.1367*irf_factors.imp_y_fac_mat_scl(7,3, :)+0.1361*irf_factors.imp_y_fac_mat_scl(8,3, :)-(irf_factors.imp_y_fac_mat_scl(6,3, :)+0.9140*irf_factors.imp_y_fac_mat_scl(7,3, :)+0.0810*irf_factors.imp_y_fac_mat_scl(8,3, :)))+0.25*squeeze(se_irf_factors.se_imp_y_fac_mat_scl(7,3, :))*2);
hold on
plot(0.25*squeeze(irf_factors.imp_y_fac_mat_scl(6,3, :)+0.1367*irf_factors.imp_y_fac_mat_scl(7,3, :)+0.1361*irf_factors.imp_y_fac_mat_scl(8,3, :)-(irf_factors.imp_y_fac_mat_scl(6,3, :)+0.9140*irf_factors.imp_y_fac_mat_scl(7,3, :)+0.0810*irf_factors.imp_y_fac_mat_scl(8,3, :)))-0.25*squeeze(se_irf_factors.se_imp_y_fac_mat_scl(7,3, :))*2);
title('Impulse Response Function: Slope (10y-3m)')
xlabel('Time period')
ylabel('Effect')
legend('Point Estimate','Upper Bound','Lower Bound')
yline(0,'--black','HandleVisibility','off')


yield_curve_irf = ns_factor_loadings(:,1)*squeeze(irf_factors.imp_y_fac_mat_scl(6,3, :))'+ns_factor_loadings(:,2)*squeeze(irf_factors.imp_y_fac_mat_scl(7,3, :))'+ns_factor_loadings(:,3)*squeeze(irf_factors.imp_y_fac_mat_scl(8,3, :))';

% Extract maturity names and time indices
maturityNames = yield_data.Properties.VariableNames;
timeIndices = 1:3:60;

% Plot the yield residuals using surf
surf(yield_curve_irf');

% Set x-axis labels to be the maturity names
xticks(1:length(maturityNames));
xticklabels(maturityNames);

yticks(timeIndices);
% Reverse the y-axis direction to go from past to recent
set(gca, 'YDir', 'reverse');

title('Yield Curve');
xlabel('Maturity');
ylabel('Time Index');
zlabel('Yields'); % Set the viewing angle for better visualization


plot(0.25*(ns_factor_loadings(1,1)*squeeze(irf_factors.imp_y_fac_mat_scl(6,3, :))'+ns_factor_loadings(1,2)*squeeze(irf_factors.imp_y_fac_mat_scl(7,3, :))'+ns_factor_loadings(1,3)*squeeze(irf_factors.imp_y_fac_mat_scl(8,3, :))'));
hold on
plot(0.25*squeeze(ns_factor_loadings(1,1)*squeeze(irf_factors.imp_y_fac_mat_scl(6,3, :))'+ns_factor_loadings(1,2)*squeeze(irf_factors.imp_y_fac_mat_scl(7,3, :))'+ns_factor_loadings(1,3)*squeeze(irf_factors.imp_y_fac_mat_scl(8,3, :))')+0.25*(ns_factor_loadings(1,1)*squeeze(se_irf_factors.se_imp_y_fac_mat_scl(6,3, :))'+ns_factor_loadings(1,2)*squeeze(se_irf_factors.se_imp_y_fac_mat_scl(7,3, :))'+ns_factor_loadings(1,3)*squeeze(se_irf_factors.se_imp_y_fac_mat_scl(8,3, :))')*2);
hold on
plot(0.25*squeeze(ns_factor_loadings(1,1)*squeeze(irf_factors.imp_y_fac_mat_scl(6,3, :))'+ns_factor_loadings(1,2)*squeeze(irf_factors.imp_y_fac_mat_scl(7,3, :))'+ns_factor_loadings(1,3)*squeeze(irf_factors.imp_y_fac_mat_scl(8,3, :))')-0.25*(ns_factor_loadings(1,1)*squeeze(se_irf_factors.se_imp_y_fac_mat_scl(6,3, :))'+ns_factor_loadings(1,2)*squeeze(se_irf_factors.se_imp_y_fac_mat_scl(7,3, :))'+ns_factor_loadings(1,3)*squeeze(se_irf_factors.se_imp_y_fac_mat_scl(8,3, :))')*2);
title('Impulse Response Function: 3-month yield')
xlabel('Time period')
ylabel('Effect')
legend('Point Estimate','Upper Bound','Lower Bound')
yline(0,'--black','HandleVisibility','off')


plot(0.25*(ns_factor_loadings(4,1)*squeeze(irf_factors.imp_y_fac_mat_scl(6,3, :))'+ns_factor_loadings(4,2)*squeeze(irf_factors.imp_y_fac_mat_scl(7,3, :))'+ns_factor_loadings(4,3)*squeeze(irf_factors.imp_y_fac_mat_scl(8,3, :))'));
hold on
plot(0.25*squeeze(ns_factor_loadings(4,1)*squeeze(irf_factors.imp_y_fac_mat_scl(6,3, :))'+ns_factor_loadings(4,2)*squeeze(irf_factors.imp_y_fac_mat_scl(7,3, :))'+ns_factor_loadings(4,3)*squeeze(irf_factors.imp_y_fac_mat_scl(8,3, :))')+0.25*(ns_factor_loadings(4,1)*squeeze(se_irf_factors.se_imp_y_fac_mat_scl(6,3, :))'+ns_factor_loadings(4,2)*squeeze(se_irf_factors.se_imp_y_fac_mat_scl(7,3, :))'+ns_factor_loadings(4,3)*squeeze(se_irf_factors.se_imp_y_fac_mat_scl(8,3, :))')*2);
hold on
plot(0.25*squeeze(ns_factor_loadings(4,1)*squeeze(irf_factors.imp_y_fac_mat_scl(6,3, :))'+ns_factor_loadings(4,2)*squeeze(irf_factors.imp_y_fac_mat_scl(7,3, :))'+ns_factor_loadings(4,3)*squeeze(irf_factors.imp_y_fac_mat_scl(8,3, :))')-0.25*(ns_factor_loadings(4,1)*squeeze(se_irf_factors.se_imp_y_fac_mat_scl(6,3, :))'+ns_factor_loadings(4,2)*squeeze(se_irf_factors.se_imp_y_fac_mat_scl(7,3, :))'+ns_factor_loadings(4,3)*squeeze(se_irf_factors.se_imp_y_fac_mat_scl(8,3, :))')*2);
title('Impulse Response Function: 2-year yield')
xlabel('Time period')
ylabel('Effect')
legend('Point Estimate','Upper Bound','Lower Bound')
yline(0,'--black','HandleVisibility','off')


plot(0.25*(ns_factor_loadings(10,1)*squeeze(irf_factors.imp_y_fac_mat_scl(6,3, :))'+ns_factor_loadings(10,2)*squeeze(irf_factors.imp_y_fac_mat_scl(7,3, :))'+ns_factor_loadings(10,3)*squeeze(irf_factors.imp_y_fac_mat_scl(8,3, :))'));
hold on
plot(0.25*squeeze(ns_factor_loadings(10,1)*squeeze(irf_factors.imp_y_fac_mat_scl(6,3, :))'+ns_factor_loadings(10,2)*squeeze(irf_factors.imp_y_fac_mat_scl(7,3, :))'+ns_factor_loadings(10,3)*squeeze(irf_factors.imp_y_fac_mat_scl(8,3, :))')+0.25*(ns_factor_loadings(10,1)*squeeze(se_irf_factors.se_imp_y_fac_mat_scl(6,3, :))'+ns_factor_loadings(10,2)*squeeze(se_irf_factors.se_imp_y_fac_mat_scl(7,3, :))'+ns_factor_loadings(10,3)*squeeze(se_irf_factors.se_imp_y_fac_mat_scl(8,3, :))')*2);
hold on
plot(0.25*squeeze(ns_factor_loadings(10,1)*squeeze(irf_factors.imp_y_fac_mat_scl(6,3, :))'+ns_factor_loadings(10,2)*squeeze(irf_factors.imp_y_fac_mat_scl(7,3, :))'+ns_factor_loadings(10,3)*squeeze(irf_factors.imp_y_fac_mat_scl(8,3, :))')-0.25*(ns_factor_loadings(10,1)*squeeze(se_irf_factors.se_imp_y_fac_mat_scl(6,3, :))'+ns_factor_loadings(10,2)*squeeze(se_irf_factors.se_imp_y_fac_mat_scl(7,3, :))'+ns_factor_loadings(10,3)*squeeze(se_irf_factors.se_imp_y_fac_mat_scl(8,3, :))')*2);
title('Impulse Response Function: 8-year yield')
xlabel('Time period')
ylabel('Effect')
legend('Point Estimate','Upper Bound','Lower Bound')
yline(0,'--black','HandleVisibility','off')

