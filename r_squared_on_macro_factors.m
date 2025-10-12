%R-squared of observables on macro factors

R_squared_on_macro_factors_slow = zeros(0, size(slow_moving_variables_final,2));
R_squared_on_macro_factors_fast = zeros(0, size(fast_moving_variables_final,2));

for i=1:size(slow_moving_variables_final,2)

Y=table2array(slow_moving_variables_final(:,i));

X=table2array([table2array(slow_macro_factors)]);

% Fit the linear regression model
mdl = fitlm(X, Y);

% Display the model summary (optional)
disp(mdl);

% Extract the R-squared value
R_squared_on_macro_factors_slow(i) = mdl.Rsquared.Ordinary;
% Convert the array to a table with column names
end


R_squared_on_macro_factors_slow = array2table(R_squared_on_macro_factors_slow, 'VariableNames', slow_moving_variables_final.Properties.VariableNames);

R_squared_on_macro_factors_slow.Properties.RowNames = "R-squared";

% Convert table to cell array including row names
R_squared_on_macro_factors_slow = [ [{'RowNames'}, R_squared_on_macro_factors_slow.Properties.VariableNames];  % Header
               [R_squared_on_macro_factors_slow.Properties.RowNames, table2cell(R_squared_on_macro_factors_slow)] ];   % Data with row names

% Export to Excel
writecell(R_squared_on_macro_factors_slow, 'R_squared_on_macro_factors_slow.xlsx');


for i=1:size(fast_moving_variables_final,2)

Y=table2array(fast_moving_variables_final(:,i));

X=table2array([table2array(slow_macro_factors) table2array(fast_macro_factors)]);

% Fit the linear regression model
mdl = fitlm(X, Y);

% Display the model summary (optional)
disp(mdl);

% Extract the R-squared value
R_squared_on_macro_factors_fast(i) = mdl.Rsquared.Ordinary;
% Convert the array to a table with column names
end

R_squared_on_macro_factors_fast = array2table(R_squared_on_macro_factors_fast, 'VariableNames', fast_moving_variables_final.Properties.VariableNames);

R_squared_on_macro_factors_fast.Properties.RowNames = "R-squared";

% Convert table to cell array including row names
R_squared_on_macro_factors_fast = [ [{'RowNames'}, R_squared_on_macro_factors_fast.Properties.VariableNames];  % Header
               [R_squared_on_macro_factors_fast.Properties.RowNames, table2cell(R_squared_on_macro_factors_fast)] ];   % Data with row names

% Export to Excel
writecell(R_squared_on_macro_factors_fast, 'R_squared_on_macro_factors_fast.xlsx');



