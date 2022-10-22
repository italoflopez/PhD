%bla
yield_data=xlsread("C:\Users\Italo\Documents\PhD\PhD\USD Zero Coupon Yields.xlsx");
macro_data=readtable("C:\Users\Italo\Documents\PhD\PhD\final_macro_data.csv");
lambda=0.0609;
maturity=[3, 6, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 180, 240, 360];
factor_loading_beta1 = ones(1,length(maturity));
factor_loading_beta2 = (1-exp(-lambda*maturity))./(lambda*maturity);
factor_loading_beta3 = (1-exp(-lambda*maturity))./(lambda*maturity)-exp(-lambda*maturity);
ns_factor_loadings = [ones(1,length(maturity)); (1-exp(-lambda*maturity))./(lambda*maturity); (1-exp(-lambda*maturity))./(lambda*maturity)-exp(-lambda*maturity)]';
for i=1:height(yield_data);
betas(i,:)=inv(ns_factor_loadings'*ns_factor_loadings)*ns_factor_loadings'*yield_data(i,:)';
end

[h,pValue,stat,cValue] = adftest(betas(:,1))
[h,pValue,stat,cValue] = adftest(betas(:,2))
[h,pValue,stat,cValue] = adftest(betas(:,3))

[h,pValue,stat,cValue] = kpsstest(betas(:,1))
[h,pValue,stat,cValue] = kpsstest(betas(:,2))
[h,pValue,stat,cValue] = kpsstest(betas(:,3))

macro_data=readtable("C:\Users\Italo\Documents\PhD\PhD\standardize_macro_data.csv");

sigma_x=cov(table2array(macro_data(:,2:44)));
[V,D] = eig(sigma_x);

big_lambda=V(:,((length(V)-4):length(V)));
big_lambda=big_lambda(:,[5 4 3 2 1]);

ft=table2array(macro_data(:,2:44))*big_lambda/height(macro_data(:,2:44));

common_component=(big_lambda*ft')';

idiosyncratic_error=table2array(macro_data(:,2:44))-common_component;

for i=1:43
macro_factor_loadings(:,i)=inv(ft'*ft)*ft'*table2array(macro_data(:,(1+i)));
end

common_component=ft*macro_factor_loadings;

idiosyncratic_error=table2array(macro_data(:,2:44))-common_component;

macro_factor_loadings*macro_factor_loadings'/height(macro_data(:,2:44));
diag(D/sum(diag(D)))*100;

macro_data = [macro_data(:,1) macro_data(:,31) macro_data(:,2:30) macro_data(:,32:50)];

macro_data = [macro_data(:,1:39) macro_data(:,41) macro_data(:,40) macro_data(:,42:50)];

macro_data = [macro_data(:,1:40) macro_data(:,48) macro_data(:,41) macro_data(:,42:47) macro_data(:,49:50)];

macro_data = [macro_data(:,1:4) macro_data(:,41) macro_data(:,5:40) macro_data(:,42:50)];

macro_data = [macro_data(:,1:41) macro_data(:,49) macro_data(:,42:48) macro_data(:,50)];

slow_moving_variables=table2array(macro_data(:,2:40));

fast_moving_variables=table2array(macro_data(:,41:50));

sigma_x_slow=cov(slow_moving_variables);

sigma_x_fast=cov(fast_moving_variables);

[V_slow,D_slow] = eig(sigma_x_slow);

[V_fast,D_fast] = eig(sigma_x_fast);

big_lambda_slow=V_slow(:,((length(V_slow)-1):length(V_slow)));
big_lambda_slow=big_lambda_slow(:,[2 1]);

big_lambda_fast=V_fast(:,((length(V_fast)-1):length(V_fast)));
big_lambda_fast=big_lambda_fast(:,[2 1]);

ft_slow=slow_moving_variables*big_lambda_slow/height(slow_moving_variables);

ft_fast=fast_moving_variables*big_lambda_fast/height(fast_moving_variables);

diag(D_slow/sum(diag(D_slow)))*100;

diag(D_fast/sum(diag(D_fast)))*100;
