% my_functions_script.m
function betas=betas_function(yield_data)

lambda=0.0609;
maturity=[3, 6, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 180, 240, 360];
factor_loading_beta1 = ones(1,length(maturity));
factor_loading_beta2 = (1-exp(-lambda*maturity))./(lambda*maturity);
factor_loading_beta3 = (1-exp(-lambda*maturity))./(lambda*maturity)-exp(-lambda*maturity);
%Nelson-Siegel Factor Loadings
ns_factor_loadings = [ones(1,length(maturity)); (1-exp(-lambda*maturity))./(lambda*maturity); (1-exp(-lambda*maturity))./(lambda*maturity)-exp(-lambda*maturity)]';
%Nelson-Siegel Factors
matrix_inverted_ns=inv(ns_factor_loadings'*ns_factor_loadings);
for i=1:height(yield_data);
betas(i,:)=matrix_inverted_ns*ns_factor_loadings'*table2array(yield_data(i,:))';
end

betas=table2timetable([table(datetime(yield_data.Var1,'Format','yyyyMM')), array2table(betas)]);

function table_betas=table_betas_function(betas)
table_betas=array2table([mean(table2array(betas))' min(table2array(betas))' max(table2array(betas))' std(table2array(betas))' [adftest(betas(:,1)).pValue adftest(betas(:,2)).pValue adftest(betas(:,3)).pValue]'])
table_betas.Properties.VariableNames = ["Mean","Minimum","Maximum","Satndard Deviation","ADF test p-value"]
table_betas=[array2table(["Beta 1","Beta 2","Beta 3"]') table_betas];
table_betas.Properties.VariableNames(1)="Factor"
table_betas

