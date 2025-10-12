function table_betas=table_betas_function(betas)
table_betas=array2table([mean(table2array(betas))' min(table2array(betas))' max(table2array(betas))' std(table2array(betas))' [adftest(betas(:,1)).pValue adftest(betas(:,2)).pValue adftest(betas(:,3)).pValue]'])
table_betas.Properties.VariableNames = ["Mean","Minimum","Maximum","Satndard Deviation","ADF test p-value"]
table_betas=[array2table(["Beta 1","Beta 2","Beta 3"]') table_betas];
table_betas.Properties.VariableNames(1)="Factor"
table_betas



