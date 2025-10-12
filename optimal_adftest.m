% Example data
y = out.fac(:,8); % Replace with your time series data
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

% Display the results
disp(['Optimal lag length (based on AIC): ', num2str(optimalLag - 1)]);
disp('AIC values for each lag length:');
disp(aicValues);

% Perform ADF test with the optimal lag length
[h, pValue, stat, cValue] = adftest(y, 'Lags', optimalLag - 1);

% Display ADF results
disp(['ADF Test Hypothesis Rejected (h=1 means reject): ', num2str(h)]);
disp(['p-Value: ', num2str(pValue)]);
disp(['Test Statistic: ', num2str(stat)]);
disp(['Critical Value: ', num2str(cValue)]);