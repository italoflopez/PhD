cd("C:\Users\Italo\Documents\PhD\PhD\First")

%Retrieveing yield data
yield_data=readtable("C:\Users\Italo\Documents\PhD\PhD\USD Zero Coupon Yields.xlsx");

yield_data=yield_data(3:height(yield_data),:);


%Putting yield data in timetable format
t1 = datetime(1989,3,1);
t2 = datetime(2022,6,1);
t = t1:calmonths(1):t2;
t3 = datetime(1989,3,1);
t4 = datetime(2007,12,1);
t_l = t3:calmonths(1):t4;
t5 = datetime(2015,1,1);
t_last = t5:calmonths(1):t2;
start_t_last=t1:calmonths(1):t5;



yield_data=table2timetable([table(datetime(t_l','Format','yyyyMM')), yield_data(1:size(t_l',1),2:size(yield_data,2))]);


[betas,ns_factor_loadings] = betas_function(yield_data);

%Plots of the betas
plot(betas.Var1,betas(:,1).betas1);
title('Beta 1')
plot(betas.Var1,betas(:,2).betas2);
title('Beta 2')
plot(betas.Var1,betas(:,3).betas3);
title('Beta 3')

%stationarity tests
[h,pValue] = adftest(betas(:,1));
[h,pValue] = adftest(betas(:,2));
[h,pValue] = adftest(betas(:,3));

[h,pValue] = kpsstest(betas(:,1));
[h,pValue] = kpsstest(betas(:,2));
[h,pValue] = kpsstest(betas(:,3));

[h,pValue] = pptest(betas(:,1));
[h,pValue] = pptest(betas(:,2));
[h,pValue] = pptest(betas(:,3));

table_betas=table_betas_function(betas)
writetable(table_betas,'betas_table.xlsx');


stationary_betas=array2table([diff(table2array(betas(:,1:2))) table2array(betas(2:end,end))]);
stationary_betas = renamevars(stationary_betas,["Var1","Var2","Var3"],["Beta1","Beta2","Beta3"]); 
stationary_betas=table2timetable([table(datetime(yield_data(2:end,:).Var1,'Format','yyyyMM')), stationary_betas]);


[estimated_yield_curve,yield_residuals]=yield_residuals_function(yield_data,ns_factor_loadings,betas);
surf(table2array(yield_residuals));
title('Yield Curve Residuals');

macro_data=readtable("C:\Users\Italo\Documents\PhD\PhD\First\stationary_data_for_macro_factors.csv");
macro_data=table2timetable([table(datetime(macro_data.Var1,'Format','yyyy-MM-dd')), macro_data(:,2:size(macro_data,2))]);
[coeff,score,latent,tsquared,explained,mu] = pca(zscore(table2array(macro_data)));
pca_macro_factors=zscore(table2array(macro_data))*coeff;
pca_macro_factors=array2timetable(pca_macro_factors,'RowTimes', datetime(macro_data.Var1,'Format','yyyyMM'));


plot(pca_macro_factors.Time,pca_macro_factors.pca_macro_factors1);
hold on 
plot(pca_macro_factors.Time,pca_macro_factors.pca_macro_factors2);
hold on
plot(pca_macro_factors.Time,pca_macro_factors.pca_macro_factors3);
hold on
plot(pca_macro_factors.Time,pca_macro_factors.pca_macro_factors4);

proportion_of_variances = std(table2array(pca_macro_factors))/sum(std(table2array(pca_macro_factors)))*100
sum(proportion_of_variances)
cumsum(proportion_of_variances)
sum(proportion_of_variances(1:4))
pca_macro_factors=pca_macro_factors(:,1:4);


ffr=readtable("C:\Users\Italo\Documents\PhD\PhD\ffr.csv");
ffr=table2timetable([table(datetime(char(ffr{:,1}),'Format','yyyyMM')), ffr(:,2:size(ffr,2))]);

diff_ffr = diff(table2array(ffr));
diff_ffr = array2timetable(diff_ffr,'RowTimes',ffr(2:end,:).Var1);

data=synchronize(macro_data,ffr,betas);
data=rmmissing(data);
macro_data=data(:,1:78);

slow_moving_variables=macro_data(:,2:47);
slow_moving_variables=[slow_moving_variables macro_data(:,59)];
slow_moving_variables=[slow_moving_variables macro_data(:,61:78)];


[coeff,score,latent,tsquared,explained,mu] = pca(zscore(table2array(slow_moving_variables)));
slow_macro_factors=zscore(table2array(slow_moving_variables))*coeff;
number_of_slow_factors=2;
slow_macro_factors=slow_macro_factors(:,1:number_of_slow_factors);
slow_macro_factors=array2timetable(array2table(slow_macro_factors),'RowTimes', datetime(macro_data.Var1,'Format','yyyyMM'));
std(table2array(table2array(slow_macro_factors)))


bar(explained)
ylabel('Marginal R2')
xlabel('Number of factors')
legend('Marginal R2 for slow-moving variables', 'Location', 'northeast')

plot(slow_macro_factors.Time,table2array(slow_macro_factors.Var1));
hold on 
plot(slow_macro_factors.Time,table2array(slow_macro_factors.Var2));


fast_moving_variables=macro_data(:,48:58);

[coeff,score,latent,tsquared,explained,mu] = pca(zscore(table2array(fast_moving_variables)));
fast_macro_factors=zscore(table2array(fast_moving_variables))*coeff;
number_of_fast_macro_factors=2;
fast_macro_factors=fast_macro_factors(:,1:number_of_fast_macro_factors);
fast_macro_factors=array2timetable(array2table(fast_macro_factors),'RowTimes', datetime(macro_data.Var1,'Format','yyyyMM'));
std(table2array(table2array(fast_macro_factors)))

bar(explained)
ylabel('Marginal R2')
xlabel('Number of factors')
legend('Marginal R2 for fast-moving variables', 'Location', 'northeast')

plot(fast_macro_factors.Time,table2array(fast_macro_factors.Var1));
hold on 
plot(fast_macro_factors.Time,table2array(fast_macro_factors.Var2));

%let's reorder slow moving variables
slow_moving_variables=[slow_moving_variables(:,1) slow_moving_variables(:,48) slow_moving_variables(:,2:47) slow_moving_variables(:,49:65)];
%reordering fast moving variables
fast_moving_variables=[fast_moving_variables(:,1) fast_moving_variables(:,6) fast_moving_variables(:,2:5) fast_moving_variables(:,7:11)];
matrix_inverted_slow=inv(table2array(table2array(slow_macro_factors))'*table2array(table2array(slow_macro_factors)));


for i=1:size(slow_moving_variables,2)
pca_slow_factor_loadings(:,i)=matrix_inverted_slow*table2array(table2array(slow_macro_factors))'*table2array(slow_moving_variables(:,i));
end

matrix_inverted_fast=inv(table2array(table2array(fast_macro_factors))'*table2array(table2array(fast_macro_factors)));
for i=1:size(fast_moving_variables,2)
pca_fast_factor_loadings(:,i)=matrix_inverted_fast*table2array(table2array(fast_macro_factors))'*table2array(fast_moving_variables(:,i));
end



number_of_named_slow_macro_factors=2;
named_slow_macro_factors=table2array(table2array(slow_macro_factors))*pca_slow_factor_loadings(:,1:number_of_named_slow_macro_factors);
%corr(named_slow_macro_factors(:,1),named_slow_macro_factors(:,2));
named_slow_macro_factors=table2timetable([table(datetime(slow_moving_variables.Var1,'Format','yyyyMM')), array2table(named_slow_macro_factors)]);

number_of_named_fast_macro_factors=2;
named_fast_macro_factors=table2array(table2array(fast_macro_factors))*pca_fast_factor_loadings(:,1:number_of_named_fast_macro_factors);
%corr(named_fast_macro_factors(:,1),named_fast_macro_factors(:,2));
named_fast_macro_factors=table2timetable([table(datetime(slow_moving_variables.Var1,'Format','yyyyMM')), array2table(named_fast_macro_factors)]);

regression_factors=synchronize(named_slow_macro_factors,diff_ffr,named_fast_macro_factors);
regression_factors=rmmissing(regression_factors);

full_sample=synchronize(named_slow_macro_factors,diff_ffr,named_fast_macro_factors,slow_moving_variables,fast_moving_variables);
full_sample=rmmissing(full_sample);
named_slow_macro_factors=full_sample(:,1:2);
slow_moving_variables=full_sample(:,6:70);
constant=ones(size(named_slow_macro_factors,1),1);
matrix_inverted_named_slow=inv([table2array(named_slow_macro_factors) constant]'*[table2array(named_slow_macro_factors) constant]);
for i=1:size(slow_moving_variables,2)
named_slow_factor_loadings(:,i)=matrix_inverted_named_slow*[table2array(named_slow_macro_factors) constant]'*[table2array(slow_moving_variables(:,i))];
end

for i=1:size(slow_moving_variables,2)
slow_variables_hat(:,i)=[table2array(named_slow_macro_factors) constant]*named_slow_factor_loadings(:,i);
end
slow_variables_residuals=table2array(slow_moving_variables)-slow_variables_hat;
slow_variables_residuals=array2timetable(slow_variables_residuals,'RowTimes',datetime(slow_moving_variables.Var1,'Format','yyyyMM'));

constant=ones(size(regression_factors,1),1);
fast_moving_variables=full_sample(:,71:81);
matrix_inverted_named_fast = inv([table2array(regression_factors) constant]'*[table2array(regression_factors) constant]);
for i=1:size(fast_moving_variables,2)
named_fast_factor_loadings(:,i)=matrix_inverted_named_fast*[table2array(regression_factors) constant]'*table2array(fast_moving_variables(:,i));
end

for i=1:size(fast_moving_variables,2)
fast_variables_hat(:,i)=[table2array(regression_factors) constant]*named_fast_factor_loadings(:,i);
end
fast_variables_residuals=table2array(fast_moving_variables)-fast_variables_hat;
fast_variables_residuals=array2timetable(fast_variables_residuals,'RowTimes',datetime(fast_moving_variables.Var1,'Format','yyyyMM'));

factors=synchronize(regression_factors,stationary_betas);
factors=rmmissing(factors);

number_of_factors = size(factors,2);
lag_order = 5;
Mdl = varm(number_of_factors,lag_order);
%EstMdl = estimate(Mdl,table2array(factors),X=double([dummy1 dummy2]));
% factors = filloutliers(factors,"linear");
EstMdl = estimate(Mdl,factors);
var_summary=summarize(EstMdl)
resid_covariance=var_summary.Covariance;
H=chol(resid_covariance,'lower')*diag(diag(chol(resid_covariance,'lower')).^2)^(-1/2);
E = infer(EstMdl,factors);
response=irf(EstMdl);
armairf(EstMdl.AR,[],InnovCov=var_summary.Covariance);
%Estimation of the IRFs
A_companion_form =[EstMdl.AR{1,1:lag_order};eye(number_of_factors*(lag_order-1)) zeros(number_of_factors*(lag_order-1),number_of_factors)];% 
eig_companion_form=eig(A_companion_form);
abs(eig_companion_form);
AICs=lag_order_selection(factors,12);
[min_num,min_idx] = min(AICs);
for i = 1:number_of_factors
table_residuals_autocorrelation_test(i,:) = lbqtest(E(:,(i+number_of_factors)));
end
table_residuals_autocorrelation_test
[h,pValue,stat,cValue]=mlbqtest(table2array(E(:,(number_of_factors+1):(2*number_of_factors))),24);

values = ffr.EffectiveFederalFundsRate;

belowThreshold = ffr(values < 0.25, :);