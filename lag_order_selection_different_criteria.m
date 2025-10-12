% Define maximum lag
maxLag = 6;  

%variables

Y = out.fac;

[n, numVars] = size(Y);

% Initialize storage for AIC, BIC, HQIC
aic = zeros(maxLag,1);
bic = zeros(maxLag,1);
hqic = zeros(maxLag,1);

for p = 1:maxLag
    % Estimate VAR(p) model
    Mdl = varm(numVars, p);
    EstMdl = estimate(Mdl, Y);
    var=summarize(EstMdl);
    
    % Extract log-likelihood
    logL = var.LogLikelihood;
    
    % Number of estimated parameters
    numParams = numel(var.Table) + numel(var.Covariance);
    
    % Compute criteria
    aic(p) = -2*logL + 2*numParams;
    bic(p) = -2*logL + numParams*log(n);
    hqic(p) = -2*logL + 2*numParams*log(log(n));
end

% Get optimal lags
[~, optAIC] = min(aic);
[~, optBIC] = min(bic);
[~, optHQIC] = min(hqic);

% Display results
disp(['Optimal lag (AIC): ', num2str(optAIC)]);
disp(['Optimal lag (BIC): ', num2str(optBIC)]);
disp(['Optimal lag (HQIC): ', num2str(optHQIC)]);
