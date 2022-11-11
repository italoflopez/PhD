%bla
yield_data=readtable("C:\Users\Italo\Documents\PhD\PhD\USD Zero Coupon Yields.xlsx");
yield_data=yield_data(3:height(yield_data),:);
%yield_data=table2timetable(yield_data);
t1 = datetime(1989,3,1);
t2 = datetime(2022,6,1);
t = t1:calmonths(1):t2
yield_data=table2timetable([table(datetime(t','Format','yyyyMM')), yield_data(:,2:size(yield_data,2))]);
%macro_data=readtable("C:\Users\Italo\Documents\PhD\PhD\final_macro_data.csv");
lambda=0.0609;
maturity=[3, 6, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 180, 240, 360];
factor_loading_beta1 = ones(1,length(maturity));
factor_loading_beta2 = (1-exp(-lambda*maturity))./(lambda*maturity);
factor_loading_beta3 = (1-exp(-lambda*maturity))./(lambda*maturity)-exp(-lambda*maturity);
%Nelson-Siegel Factor Loadings
ns_factor_loadings = [ones(1,length(maturity)); (1-exp(-lambda*maturity))./(lambda*maturity); (1-exp(-lambda*maturity))./(lambda*maturity)-exp(-lambda*maturity)]';
%Nelson-Siegel Factors
for i=1:height(yield_data);
betas(i,:)=inv(ns_factor_loadings'*ns_factor_loadings)*ns_factor_loadings'*table2array(yield_data(i,:))';
end
betas=table2timetable([table(datetime(yield_data.Var1,'Format','yyyyMM')), array2table(betas)]);
%stationarity tests
[h,pValue,stat,cValue] = adftest(betas(:,1));
[h,pValue,stat,cValue] = adftest(betas(:,2));
[h,pValue,stat,cValue] = adftest(betas(:,3));

[h,pValue,stat,cValue] = kpsstest(betas(:,1));
[h,pValue,stat,cValue] = kpsstest(betas(:,2));
[h,pValue,stat,cValue] = kpsstest(betas(:,3));

macro_data=readtable("C:\Users\Italo\Documents\PhD\PhD\First\standardize_macro_data.csv");
macro_data=table2timetable([table(datetime(macro_data.Var1,'Format','yyyyMM')), macro_data(:,2:size(macro_data,2))]);
%Variance-covariance matrix of standardized stationary macro data
sigma_x=cov(table2array(macro_data));
%eigenvalue-eigenvector decomposition
[V,D] = eig(sigma_x);
%Eigenvectors corresponding to 5 largest eigenvalues
big_lambda=V(:,((length(V)-4):length(V)));
big_lambda=big_lambda(:,[5 4 3 2 1]);
%First Macro Factors
ft=table2array(macro_data)*big_lambda/height(macro_data);
%common component
common_component=(big_lambda*ft')';
%idyosincratic eror
idiosyncratic_error=table2array(macro_data)-common_component;
%macro factor loadings by OLS
for i=1:size(macro_data,2)
macro_factor_loadings(:,i)=inv(ft'*ft)*ft'*table2array(macro_data(:,i));
end
%common component
common_component=ft*macro_factor_loadings;
%idyosincratic error
idiosyncratic_error=table2array(macro_data)-common_component;
%some checks
macro_factor_loadings*macro_factor_loadings'/height(macro_data);
diag(D/sum(diag(D)))*100;
%reordering macro data
macro_data = [macro_data(:,1) macro_data(:,30) macro_data(:,2:29) macro_data(:,31:49)];

macro_data = [macro_data(:,1:38) macro_data(:,40) macro_data(:,39) macro_data(:,41:49)];

macro_data = [macro_data(:,1:39) macro_data(:,46) macro_data(:,40:45) macro_data(:,47:49)];

macro_data = [macro_data(:,1:2) macro_data(:,47) macro_data(:,3:46) macro_data(:,48:49)];

%macro_data = [macro_data(:,1:41) macro_data(:,49) macro_data(:,42:48)];
%slow moving variables
slow_moving_variables=macro_data(:,1:39);
%fast moving variables
fast_moving_variables=macro_data(:,40:49);
%variance-covariance matrix of slow moving variables
sigma_x_slow=cov(table2array(slow_moving_variables));
%variance-covariance matrix of fast moving variables
sigma_x_fast=cov(table2array(fast_moving_variables));
%eingenvalue-eigenvector decomposition
[V_slow,D_slow] = eig(sigma_x_slow);
%eigenvalue-eigenvector decomposition
[V_fast,D_fast] = eig(sigma_x_fast);

big_lambda_slow=V_slow(:,((length(V_slow)-1):length(V_slow)));
big_lambda_slow=big_lambda_slow(:,[2 1]);

big_lambda_fast=V_fast(:,((length(V_fast)-1):length(V_fast)));
big_lambda_fast=big_lambda_fast(:,[2 1]);
%%%%%%%
ft_slow=table2array(slow_moving_variables)*big_lambda_slow/height(slow_moving_variables);
ft_slow=table2timetable([table(datetime(slow_moving_variables.Var1,'Format','yyyyMM')), array2table(ft_slow)]);
ft_fast=table2array(fast_moving_variables)*big_lambda_fast/height(fast_moving_variables);
ft_fast=table2timetable([table(datetime(fast_moving_variables.Var1,'Format','yyyyMM')), array2table(ft_fast)]);

diag(D_slow/sum(diag(D_slow)))*100;

diag(D_fast/sum(diag(D_fast)))*100;
%factor loadings of slow factors on slow moving variables
for i=1:size(slow_moving_variables,2)
slow_factor_loadings(:,i)=inv(table2rray(ft_slow)'*table2array(ft_slow))*table2array(ft_slow)'*table2array(slow_moving_variables(:,i));
end
%factor loadings of fast factors on macro variables
for i=1:size(macro_data,2)
fast_factor_loadings(:,i)=inv(table2array(ft_fast)'*table2array(ft_fast))*table2array(ft_fast)'*table2array(macro_data(:,i));
end
%slow factors with name normalization
ft_slow_name=table2array(ft_slow)*slow_factor_loadings(:,1:2);
ft_slow_name=table2timetable([table(datetime(slow_moving_variables.Var1,'Format','yyyyMM')), array2table(ft_slow_name)]);
%fast factors with name normalization
ft_fast_name=table2array(ft_fast)*fast_factor_loadings(:,40:41);
ft_fast_name=table2timetable([table(datetime(fast_moving_variables.Var1,'Format','yyyyMM')), array2table(ft_fast_name)]);

ffr=readtable("C:\Users\Italo\Documents\PhD\PhD\ffr.csv");
ffr=table2timetable([table(datetime(char(ffr{:,1}),'Format','yyyyMM')), ffr(:,2:size(ffr,2))]);

factors=synchronize(ft_slow_name,ffr,ft_fast_name,betas);
factors=rmmissing(factors);

Mdl = varm(8,2);
EstMdl = estimate(Mdl,factors);
var=summarize(EstMdl)
resid_covariance=var.Covariance;
H=chol(resid_covariance,'lower')*diag(diag(chol(resid_covariance,'lower')).^2)^(-1/2);
E = infer(EstMdl,factors);

fast_and_slow_variables=[slow_moving_variables fast_moving_variables];
r=zeros(size(fast_and_slow_variables,2),1);

b=ones(size(slow_moving_variables,2),1);
a=zeros(size(slow_moving_variables,2),2);
c=ones(size(slow_moving_variables,2),5);
f=zeros(size(fast_moving_variables,2),8);
R=[a b c;f];

for i=1:size(fast_and_slow_variables,2)
observation_equation_factor_loadings(:,i)=inv(table2array(factors)'*table2array(factors))*table2array(factors)'*table2array(fast_and_slow_variables(:,i));
end

P=R*inv(table2array(factors)'*table2array(factors))*R';
factor_loadings_with_restrictions=observation_equation_factor_loadings-inv(table2array(factors)'*table2array(factors))*R'*inv(P)*(R*observation_equation_factor_loadings-r);

small_R=[0 0 1 0 0 0 0 0;
         0 0 0 1 0 0 0 0;
         0 0 0 0 1 0 0 0;
         0 0 0 0 0 1 0 0;
         0 0 0 0 0 0 1 0;
         0 0 0 0 0 0 0 1];
small_r=[0;0;0;0;0;0];
small_P=small_R*inv(table2array(factors)'*table2array(factors))*small_R';

for i=1:size(slow_moving_variables,2)
small_factor_loadings_with_restrictions(:,i)=observation_equation_factor_loadings(:,i)-inv(table2array(factors)'*table2array(factors))*small_R'*inv(small_P)*(small_R*observation_equation_factor_loadings(:,i)-small_r);
end
