cd("C:\Users\Italo\OneDrive\Documents\PhD\PhD\First")
addpath('C:\Users\Italo\OneDrive\Documents\PhD\PhD\First\ddisk\matlab');
%Retrieveing yield data
yield_data=readtable("C:\Users\Italo\OneDrive\Documents\PhD\PhD\USD Zero Coupon Yields.xlsx");

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



yield_data=table2timetable([table(datetime(t','Format','yyyyMM')), yield_data(1:size(t',1),2:size(yield_data,2))]);


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

icm=readtable("C:\Users\Italo\Documents\PhD\PhD\First\Shadow_FFR_1221.xlsx");
icm=table2timetable([table(datetime(char(icm{:,1}),'Format','yyyyMM')), icm(:,2:size(icm,2))]);

values = ffr.EffectiveFederalFundsRate;

cutoffDate = datetime(icm.Var1(end));

belowThreshold = ffr(values < 0.25 & ffr.Var1 < cutoffDate, :);

ffr_icm = ffr;

for_dates = ffr_icm(values < 0.25 & ffr_icm.Var1 < cutoffDate, :)

valuesForDesiredDate = icm.Shadow_FFR(ismember(icm.Var1, for_dates.Var1));

ffr_icm(values < 0.25 & ffr_icm.Var1 < cutoffDate, :) = array2table(valuesForDesiredDate);

data=synchronize(macro_data,ffr_icm,betas);
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

regression_factors=synchronize(named_slow_macro_factors,ffr_icm,named_fast_macro_factors);
regression_factors=rmmissing(regression_factors);

full_sample=synchronize(named_slow_macro_factors,ffr_icm,named_fast_macro_factors,slow_moving_variables,fast_moving_variables);
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

factors=synchronize(regression_factors,betas);
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


for i = 1:size(factors,2)
   plot(table2array(factors(:,i)))
   hold on
end
C(1:number_of_factors,1:number_of_factors)=eye(number_of_factors);
years=5;
months=years*12;
for i=1:(months-1)
    bla=A_companion_form^i;
C(1:number_of_factors,(number_of_factors*i+1):(number_of_factors*(i+1)))=bla(1:number_of_factors,1:number_of_factors);
end

inv_H=inv(H);
B(1:number_of_factors,1:number_of_factors)=inv_H;

for i=1:(months-1)
B(1:number_of_factors,(number_of_factors*i+1):(number_of_factors*(i+1)))=C(1:number_of_factors,(number_of_factors*i+1):(number_of_factors*(i+1)))*inv_H;
end

n_yields = size(yield_data,2);
n_non_yield_factors = size(factors,2)-size(betas,2);
n_slow_moving_variables = size(slow_moving_variables,2);
n_fast_moving_variables = size(fast_moving_variables,2);

yields_loadings=[zeros(15,5), ns_factor_loadings];
full_named_slow_factor_loadings=[zeros(5,65); named_slow_factor_loadings];
full_named_fast_factor_loadings=[zeros(2,11); named_fast_factor_loadings];
complete_observation_equation_factor_loadings=[full_named_slow_factor_loadings full_named_fast_factor_loadings, yields_loadings'];
observed_variables_irf=complete_observation_equation_factor_loadings'*B;

plot(B(3,3:number_of_factors:480));
plot(B(6,3:number_of_factors:480));
plot(B(7,3:number_of_factors:480));
plot(B(8,3:number_of_factors:480));
plot(observed_variables_irf(79,3:8:480));
%Now we will do Yamamoto (2016)
factors_starting_values=factors(1:lag_order,:);
state_equation_residuals=E(:,(number_of_factors+1):(number_of_factors+8));
state_equation_residuals=table2array(state_equation_residuals)-mean(table2array(state_equation_residuals));
n_bootstrap_samples=1000;
bootstrap_factors=zeros(size(factors,1),size(factors,2),n_bootstrap_samples);
for z=1:n_bootstrap_samples
bootstrap_sample_points=randsample((size(state_equation_residuals,1)),(size(state_equation_residuals,1)));
state_equation_bootstrap_residuals=state_equation_residuals(bootstrap_sample_points,:);
bootstrap_factors(1:lag_order,:,z)=table2array(factors_starting_values);
estimated_AR=[EstMdl.AR{1,1:lag_order}];
for i=1:(size(bootstrap_factors,1)-lag_order)
    internal_sum=estimated_AR(1:number_of_factors,1:number_of_factors)*bootstrap_factors((i+(lag_order-1)),:,z)';
    for j=1:(lag_order-1)
        internal_sum=internal_sum+estimated_AR(1:number_of_factors,(1+number_of_factors*j):(number_of_factors+number_of_factors*j))*bootstrap_factors((i+(lag_order-1)),:,z)';
    end
   bootstrap_factors(i+lag_order,:,z)= internal_sum+state_equation_bootstrap_residuals(i,:)';
end
end
plot(bootstrap_factors(:,6,100));
hold on
plot(bootstrap_factors(:,6,1));

observable_variables_residuals=synchronize(slow_variables_residuals,fast_variables_residuals,yield_residuals);
observable_variables_residuals=rmmissing(observable_variables_residuals);

observable_variables_residuals=table2array(observable_variables_residuals)-mean(table2array(observable_variables_residuals));
yields_loadings=[zeros(15,6), ns_factor_loadings];
full_named_slow_factor_loadings=[zeros(6,65); named_slow_factor_loadings];
full_named_fast_factor_loadings=[zeros(3,11); named_fast_factor_loadings];
complete_observation_equation_factor_loadings=[full_named_slow_factor_loadings full_named_fast_factor_loadings, yields_loadings'];

n_bootstrap_samples=1000;
bootstrap_observables=zeros(size(observable_variables_residuals,1),size(observable_variables_residuals,2),n_bootstrap_samples);
for z=1:n_bootstrap_samples
    bootstrap_sample_points=randsample((size(observable_variables_residuals,1)),(size(observable_variables_residuals,1)));
observable_variables_bootstrap_residuals=observable_variables_residuals(bootstrap_sample_points,:);
for i=1:size(observable_variables_residuals,1)
bootstrap_observables(i,:,z)=[bootstrap_factors(i,1:2,z) 1 bootstrap_factors(i,3:number_of_factors,z)]*complete_observation_equation_factor_loadings+observable_variables_bootstrap_residuals(i,:);
end
end

plot(bootstrap_observables(:,91,100));
hold on
plot(bootstrap_observables(:,91,1));

observables=synchronize(slow_moving_variables,fast_moving_variables,yield_data);
observables=rmmissing(observables);

plot(bootstrap_observables(:,50,100));
hold on
plot(bootstrap_observables(:,50,1));
hold on
plot(table2array(observables(:,50)));

tic;
for z=1:n_bootstrap_samples
Mdl = varm(number_of_factors,lag_order);
EstMdl = estimate(Mdl,bootstrap_factors(:,:,z));
var=summarize(EstMdl);
resid_covariance=var.Covariance;
H=chol(resid_covariance,'lower')*diag(diag(chol(resid_covariance,'lower')).^2)^(-1/2);
E = infer(EstMdl,bootstrap_factors(:,:,z));
A_companion_form =[EstMdl.AR{1,1:lag_order}; eye(number_of_factors*(lag_order-1)) zeros(number_of_factors*(lag_order-1),number_of_factors)];% 
eig_companion_form=eig(A_companion_form);
abs(eig_companion_form);

C(1:number_of_factors,1:number_of_factors)=eye(number_of_factors);
years=5;
months=years*12;
for i=1:(months-1)
    bla=A_companion_form^i;
C(1:number_of_factors,(number_of_factors*i+1):(number_of_factors*(i+1)))=bla(1:number_of_factors,1:number_of_factors);
end

inv_H=inv(H);
B_bootstrap(1:number_of_factors,1:number_of_factors,z)=inv_H;

for i=1:(months-1)
B_bootstrap(1:number_of_factors,(number_of_factors*i+1):(number_of_factors*(i+1)),z)=C(1:number_of_factors,(number_of_factors*i+1):(number_of_factors*(i+1)))*inv_H;
end

yields_loadings=[zeros(15,5), ns_factor_loadings];
full_named_slow_factor_loadings=[zeros(5,65); named_slow_factor_loadings];
full_named_fast_factor_loadings=[zeros(2,11); named_fast_factor_loadings];
complete_observation_equation_factor_loadings=[full_named_slow_factor_loadings full_named_fast_factor_loadings, yields_loadings'];
observed_variables_bootstrap_irf(:,:,z)=complete_observation_equation_factor_loadings'*B_bootstrap(:,:,z);

end
toc
std(B_bootstrap(3,3:8:480,:));
for i=0:59
irf_standard_deviation(i+1)=std(B_bootstrap(3,(3+8*i)+480*(1:999)));
end

plot(0.25*B(3,3+(8*(0:59))));
hold on
plot(0.25*B(3,3+(8*(0:59)))+0.25*irf_standard_deviation*2);
hold on
plot(0.25*B(3,3+(8*(0:59)))-0.25*irf_standard_deviation*2);
title('Impulse Response Function: Interest Rate')
xlabel('Time period')
ylabel('Effect')
legend('Upper Bound','Point Estimate','Lower Bound')
yline(0,'--black','HandleVisibility','off')


for i=0:59
irf_standard_deviation(i+1)=std(B_bootstrap(6,(3+8*i)+480*(1:999)));
end

plot(0.25*B(6,3+(8*(0:59))));
hold on
plot(0.25*B(6,3+(8*(0:59)))+0.25*irf_standard_deviation*2);
hold on
plot(0.25*B(6,3+(8*(0:59)))-0.25*irf_standard_deviation*2);
title('Impulse Response Function: Nelson-Siegel Yield Curve Level')
xlabel('Time period')
ylabel('Effect')
legend('Upper Bound','Point Estimate','Lower Bound')
yline(0,'--black','HandleVisibility','off')

for i=0:59
irf_standard_deviation(i+1)=std(B_bootstrap(7,(3+8*i)+480*(1:999)));
end

plot(-1*0.25*B(7,3+(8*(0:59))));
hold on
plot(-1*0.25*B(7,3+(8*(0:59)))+0.25*irf_standard_deviation*2);
hold on
plot(-1*0.25*B(7,3+(8*(0:59)))-0.25*irf_standard_deviation*2);
title('Impulse Response Function: Nelson-Siegel Yield Curve Slope')
xlabel('Time period')
ylabel('Effect')
legend('Upper Bound','Point Estimate','Lower Bound')
yline(0,'--black','HandleVisibility','off')

for i=0:59
irf_standard_deviation(i+1)=std(B_bootstrap(8,(3+8*i)+480*(1:999)));
end

plot(0.25*B(8,3+(8*(0:59))));
hold on
plot(0.25*B(8,3+(8*(0:59)))+0.25*irf_standard_deviation*2);
hold on
plot(0.25*B(8,3+(8*(0:59)))-0.25*irf_standard_deviation*2);
title('Impulse Response Function: Nelson-Siegel Yield Curve Curvature')
xlabel('Time period')
ylabel('Effect')
legend('Upper Bound','Point Estimate','Lower Bound')
yline(0,'--black','HandleVisibility','off')

for i=0:59
irf_standard_deviation(i+1)=std(observed_variables_bootstrap_irf(1,(3+8*i)+480*(1:999)));
end

plot(0.25*observed_variables_irf(1,3+(8*(0:59))));
hold on
plot(0.25*observed_variables_irf(1,3+(8*(0:59)))+0.25*irf_standard_deviation*2);
hold on
plot(0.25*observed_variables_irf(1,3+(8*(0:59)))-0.25*irf_standard_deviation*2);
title('Impulse Response Function: Industrial Production, final goods')
xlabel('Time period')
ylabel('Effect')
legend('Upper Bound','Point Estimate','Lower Bound')
yline(0,'--black','HandleVisibility','off')

for i=0:59
irf_standard_deviation(i+1)=std(observed_variables_bootstrap_irf(2,(3+8*i)+480*(1:999)));
end

plot(0.25*observed_variables_irf(2,3+(8*(0:59))));
hold on
plot(0.25*observed_variables_irf(2,3+(8*(0:59)))+0.25*irf_standard_deviation*2);
hold on
plot(0.25*observed_variables_irf(2,3+(8*(0:59)))-0.25*irf_standard_deviation*2);
title('Impulse Response Function: CPI')
xlabel('Time period')
ylabel('Effect')
legend('Upper Bound','Point Estimate','Lower Bound')
yline(0,'--black','HandleVisibility','off')

for i=0:59
irf_standard_deviation(i+1)=std(observed_variables_bootstrap_irf(26,(3+8*i)+480*(1:999)));
end

plot(0.25*observed_variables_irf(26,3+(8*(0:59))));
hold on
plot(0.25*observed_variables_irf(26,3+(8*(0:59)))+0.25*irf_standard_deviation*2);
hold on
plot(0.25*observed_variables_irf(26,3+(8*(0:59)))-0.25*irf_standard_deviation*2);
title('Impulse Response Function: Long-term Unemployment')
xlabel('Time period')
ylabel('Effect')
legend('Upper Bound','Point Estimate','Lower Bound')
yline(0,'--black','HandleVisibility','off')
