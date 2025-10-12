% variables_to_test=synchronize(macro_data_final,diffbetas);
% variables_to_test=rmmissing(variables_to_test);
% % Assume your timetable is called 'tt'
% numVars = width(variables_to_test); % Number of variables (columns)
% varNames = variables_to_test.Properties.VariableNames; % Variable names

%Without outliers
numVars = width(cleanedTimetable); % Number of variables (columns)
varNames = cleanedTimetable.Properties.VariableNames;

stationaryVars = {}; % To store names of stationary variables
nonStationaryVars = {}; % To store names of non-stationary variables

% Loop through each variable and perform the ADF test
for i = 1:numVars
    data = cleanedTimetable.(varNames{i}); % Extract the ith variable
    if any(isnan(data)) % Handle missing values
        continue; % Skip this variable if it has missing data
    end
    y = data; % Replace with your time series data
    maxLag = 10; % Maximum lags to consider
    n = length(y); % Sample size
    
    % Preallocate for AIC values
    aicValues = NaN(maxLag + 1, 1);
    
    % Loop over lag lengths
    for lag = 0:maxLag
        % Prepare lagged variables
        laggedY = lagmatrix(y, 1:lag); % Create lagged terms
        diffY = diff(y); % First difference of y
        X = [ones(length(diffY), 1), laggedY(2:end, :)]; % Constant and lagged differences
        validIdx = all(~isnan(X), 2); % Remove rows with NaNs
        
        % Fit the regression: Δy_t = β0 + β1*y_t-1 + ∑ βi*Δy_t-i + ε_t
        X = X(validIdx, :);
        yDep = diffY(validIdx);
        beta = X \ yDep; % Ordinary Least Squares (OLS)
        residuals = yDep - X * beta; % Compute residuals
        
        % Compute the AIC: AIC = n * log(SSE / n) + 2 * k
        sse = sum(residuals .^ 2); % Sum of squared residuals
        k = size(X, 2); % Number of parameters (constant + lags)
        aicValues(lag + 1) = n * log(sse / n) + 2 * k;
    end
    
    % Find the lag with the minimum AIC
    [~, optimalLag] = min(aicValues);
    % Perform the ADF test
    [h, pValue] = adftest(data, 'Lags', optimalLag - 1);%
    
    % Check if the series is stationary
    if h == 1 % h = 1 indicates stationarity
        stationaryVars{end+1} = varNames{i}; % Append to stationary list
    else
        nonStationaryVars{end+1} = varNames{i}; % Append to non-stationary list
    end
end

% Display results
disp('Stationary Variables:');
disp(stationaryVars);

disp('Non-Stationary Variables:');
disp(nonStationaryVars);

% Ensure stationaryVars is a column cell array
stationaryVars = stationaryVars(:); 

% Ensure nonStationaryVars is a column cell array
nonStationaryVars = nonStationaryVars(:); 
% Optional: Save results to a file
writetable(cell2table(stationaryVars, 'VariableNames', {'StationaryVars'}), 'StationaryVariables_intermonthly_without_outliers_zlb.csv');
writetable(cell2table(nonStationaryVars, 'VariableNames', {'NonStationaryVars'}), 'NonStationaryVariables_intermonthly_without_outliers_zlb.csv');