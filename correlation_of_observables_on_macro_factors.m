%Correlation of macro factors on observables

observable = table2array(slow_moving_variables_final);

macro_factor = table2array(table2array(slow_macro_factors));

corr_slow = corr(observable,macro_factor);

columnNames = slow_moving_variables_final.Properties.VariableNames; % Extract column names

corr_slow = array2table(corr_slow);

corr_slow.Properties.RowNames = columnNames;

% Convert table to cell array including row names
corr_slow = [ [{'RowNames'}, corr_slow.Properties.VariableNames];  % Header
               [corr_slow.Properties.RowNames, table2cell(corr_slow)] ];   % Data with row names

% Export to Excel
writecell(corr_slow, 'corr_slow.xlsx');

observable = table2array(fast_moving_variables_final);

macro_factor = table2array(table2array(fast_macro_factors));

corr_fast = corr(observable,macro_factor);

columnNames = fast_moving_variables_final.Properties.VariableNames; % Extract column names

corr_fast = array2table(corr_fast);

corr_fast.Properties.RowNames = columnNames;

% Convert table to cell array including row names
corr_fast = [ [{'RowNames'}, corr_fast.Properties.VariableNames];  % Header
               [corr_fast.Properties.RowNames, table2cell(corr_fast)] ];   % Data with row names

% Export the table to an Excel file
writecell(corr_fast, 'corr_fast.xlsx');
