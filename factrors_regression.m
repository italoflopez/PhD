partial_R_squared = zeros(0, size(macro_data_final,2));
pValues = zeros(0, size(macro_data_final,2));

for i=1:size(macro_data_final,2)

Y=table2array(macro_data_final(:,i));

X=table2array([table2array(slow_macro_factors) table2array(fast_macro_factors)]);

% Fit the linear regression model
mdl_reduced = fitlm(X, Y);

% Display the model summary (optional)
disp(mdl_reduced);

% Extract the R-squared value
R_squared_Reduced = mdl_reduced.Rsquared.Ordinary;

% Extract variables from the tables within Var1 and Var2
slow_macro_factors_var1 = slow_macro_factors.Var1{:,:};
slow_macro_factors_var2 = slow_macro_factors.Var2{:,:};
fast_macro_factors_var1 = fast_macro_factors.Var1{:,:};
fast_macro_factors_var2 = fast_macro_factors.Var2{:,:};

% Rebuild the timetables with these variables
slow_macro_factors_extracted = timetable(slow_macro_factors.Time, slow_macro_factors_var1, slow_macro_factors_var2);
fast_macro_factors_extracted = timetable(fast_macro_factors.Time, fast_macro_factors_var1, fast_macro_factors_var2);

% Synchronize the extracted timetables with betas
betas_regression_data = synchronize(slow_macro_factors_extracted, fast_macro_factors_extracted, betas);
betas_regression_data=rmmissing(betas_regression_data);
regression_data=synchronize(macro_data_final(:,i),betas_regression_data);
regression_data=rmmissing(regression_data);

Y = table2array(regression_data(:,1));
X=table2array(regression_data(:,2:end));

% Fit the linear regression model
mdl_full = fitlm(X, Y);

% Display the model summary (optional)
disp(mdl_full);

% Extract the R-squared value
R_squared_Full = mdl_full.Rsquared.Ordinary;

% Extract the residual sum of squares (RSS) from both models
RSS_full = sum((mdl_full.Residuals.Raw).^2);
RSS_reduced = sum((mdl_reduced.Residuals.Raw).^2);

% Calculate partial R-squared
partial_R_squared(i) = (RSS_reduced - RSS_full) / RSS_reduced;

% R-squared values
R2_full = mdl_full.Rsquared.Ordinary;
R2_reduced = mdl_reduced.Rsquared.Ordinary;

% Number of predictors
p_full = mdl_full.NumPredictors;
p_reduced = mdl_reduced.NumPredictors;

% Number of observations
n = mdl_full.NumObservations;

% F-statistic formula for comparing nested models
F = ((R2_full - R2_reduced) / (p_full - p_reduced)) / ((1 - R2_full) / (n - p_full - 1));

% Calculate the p-value based on the F-distribution
pValues(i) = 1 - fcdf(F, p_full - p_reduced, n - p_full - 1);
end 

partial_R_squared_table = [partial_R_squared;pValues];

partial_R_squared_table = array2table(partial_R_squared_table, 'VariableNames', macro_data_final.Properties.VariableNames);

partial_R_squared_table.Properties.RowNames = ["Partial R-squared", "p-value"];

% Convert table to cell array including row names
partial_R_squared_table = [ [{'RowNames'}, partial_R_squared_table.Properties.VariableNames];  % Header
               [partial_R_squared_table.Properties.RowNames, table2cell(partial_R_squared_table)] ];   % Data with row names

% Export to Excel
writecell(partial_R_squared_table, 'partial_R_squared_table.xlsx');

