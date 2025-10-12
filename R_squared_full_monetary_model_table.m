%%%%R-squared of Full Monetary Model

R_squared_full_monetary_model = out.r2;


% Example: Get all column names
allColumnNames = model_estimation_variables.Properties.VariableNames;

% Specify the columns to exclude
excludeColumns = {'Var1', 'betas1', 'betas2','betas3'};

% Get the column names excluding the specified ones
R_squared_full_monetary_modelNames = setdiff(allColumnNames, excludeColumns, 'stable');


R_squared_full_monetary_model_table_excel = array2table(R_squared_full_monetary_model, 'RowNames', R_squared_full_monetary_modelNames);

observable_variables = model_estimation_variables.Properties.VariableNames

toRemove = ["Var1", "betas1", "betas2", "betas3"];  % Strings to remove
% Remove elements while preserving order
observable_variables(ismember(observable_variables, toRemove)) = [];

R_squared_full_monetary_model_table_excel.Properties.RowNames = observable_variables;

% Convert table to cell array including row names
R_squared_full_monetary_model_table_excel = [ [{'RowNames'}, R_squared_full_monetary_model_table_excel.Properties.VariableNames];  % Header
               [R_squared_full_monetary_model_table_excel.Properties.RowNames, table2cell(R_squared_full_monetary_model_table_excel)] ];   % Data with row names

% Export the table to an Excel file
writecell(R_squared_full_monetary_model_table_excel, 'R_squared_full_monetary_model_table.xlsx');

