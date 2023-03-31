%First Paper code
yield_data=readtable("C:\Users\Italo\Documents\PhD\PhD\USD Zero Coupon Yields.xlsx");
yield_data=yield_data(3:height(yield_data),:);
yield_data=yield_data(1:166,:);
%yield_data=table2timetable(yield_data);
t1 = datetime(1989,3,1);
t2 = datetime(2002,12,1);
t = t1:calmonths(1):t2;
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
plot(betas.Var1,betas(:,1).betas1);
title('Beta 1')
plot(betas.Var1,betas(:,2).betas2);
title('Beta 2')
plot(betas.Var1,betas(:,3).betas3);
title('Beta 3')

[h,pValue,stat,cValue] = adftest(betas(:,1));
[h,pValue,stat,cValue] = adftest(betas(:,2));
[h,pValue,stat,cValue] = adftest(betas(:,3));

[h,pValue,stat,cValue] = kpsstest(betas(:,1));
[h,pValue,stat,cValue] = kpsstest(betas(:,2));
[h,pValue,stat,cValue] = kpsstest(betas(:,3));

[h,pValue,stat,cValue] = pptest(betas(:,1));
[h,pValue,stat,cValue] = pptest(betas(:,2));
[h,pValue,stat,cValue] = pptest(betas(:,3));

table_betas=array2table([mean(table2array(betas))' min(table2array(betas))' max(table2array(betas))' std(table2array(betas))' [adftest(betas(:,1)).pValue adftest(betas(:,2)).pValue adftest(betas(:,3)).pValue]'])
table_betas.Properties.VariableNames = ["Mean","Minimum","Maximum","Satndard Deviation","ADF test p-value"]
table_betas=[array2table(["Beta 1","Beta 2","Beta 3"]') table_betas];
table_betas.Properties.VariableNames(1)="Factor"
writetable(table_betas,'betas_table.xlsx');

estimated_yield_curve=ns_factor_loadings*table2array(betas)';
estimated_yield_curve=estimated_yield_curve';
yield_residuals=table2array(yield_data)-estimated_yield_curve;
yield_residuals=table2timetable([table(datetime(yield_data.Var1,'Format','yyyyMM')), array2table(yield_residuals)]);
yield_residuals.Properties.VariableNames=yield_data.Properties.VariableNames
surf(table2array(yield_residuals));
title('Yield Curve Residuals');
estimated_yield_curve=array2timetable(array2table(estimated_yield_curve),'RowTimes', datetime(yield_data.Var1,'Format','yyyyMM'));

macro_data=readtable("C:\Users\Italo\Documents\PhD\PhD\First\stationary_data_for_macro_factors.csv");
macro_data=table2timetable([table(datetime(macro_data.Var1,'Format','yyyy-MM-dd')), macro_data(:,2:size(macro_data,2))]);
macro_data=macro_data(1:131,:);
[coeff,score,latent,tsquared,explained,mu] = pca(zscore(table2array(macro_data)));
macro_factors=zscore(table2array(macro_data))*coeff;
macro_factors=array2timetable(array2table(macro_factors),'RowTimes', datetime(macro_data.Var1,'Format','yyyyMM'));
plot(macro_factors.Time,table2array(macro_factors.Var1));
plot(macro_factors.Time,table2array(macro_factors.Var1));
hold on 
plot(macro_factors.Time,table2array(macro_factors.Var2));
hold on
plot(macro_factors.Time,table2array(macro_factors.Var3));
hold on
plot(macro_factors.Time,table2array(macro_factors.Var4));

std(table2array(table2array(macro_factors)));
macro_factors=macro_factors(:,1:4);
%fmincon with symbolic math toolbox tryout
x=sym('x',[10 1])
matlabFunction(x.^2'*x,'vars',{x},'file','objfunction');
fmincon(@objfunction,zeros(10,1));

%let's try ols with fmincon and symbolic math toolbox
x=sym('x',[2 1]);
matlabFunction((table2array(yield_data(:,1))-[table2array(yield_data(:,2)) ones(400,1)]*x)'*(table2array(yield_data(:,1))-[table2array(yield_data(:,2)) ones(400,1)]*x),'vars',{x},'file','objfunction');
fmincon(@objfunction,zeros(2,1));

%let's try pca with fmincon and symbolic math toolbox for the macro data
% with four factors
% parpool('Threads')
% x=sym('x',[(size(macro_data,2)+size(macro_data,1)) 4]);
% tic;
% matlabFunction((table2array(macro_data)'-x(1:size(macro_data,2),1:4)*x((size(macro_data,2)+1):(end),1:4)')'*(table2array(macro_data)'-x(1:size(macro_data,2),1:4)*x((size(macro_data,2)+1):(end),1:4)'),'vars',{x},'file','objfunction');
% toc
%selecting slow moving variables
ffr=readtable("C:\Users\Italo\Documents\PhD\PhD\ffr.csv");
ffr=table2timetable([table(datetime(char(ffr{:,1}),'Format','yyyyMM')), ffr(:,2:size(ffr,2))]);
ffr=ffr(1:582,:);
data=synchronize(macro_data,ffr,betas);
data=rmmissing(data);
macro_data=data(:,1:78);

slow_moving_variables=macro_data(:,2:47);
slow_moving_variables=[slow_moving_variables macro_data(:,59)];
slow_moving_variables=[slow_moving_variables macro_data(:,61:78)];
[coeff,score,latent,tsquared,explained,mu] = pca(zscore(table2array(slow_moving_variables)));
slow_macro_factors=zscore(table2array(slow_moving_variables))*coeff;
slow_macro_factors=slow_macro_factors(:,1:2);
slow_macro_factors=array2timetable(array2table(slow_macro_factors),'RowTimes', datetime(macro_data.Var1,'Format','yyyyMM'));
std(table2array(table2array(slow_macro_factors)));

plot(slow_macro_factors.Time,table2array(slow_macro_factors.Var1));
hold on 
plot(slow_macro_factors.Time,table2array(slow_macro_factors.Var2));

fast_moving_variables=macro_data(:,48:58);

[coeff,score,latent,tsquared,explained,mu] = pca(zscore(table2array(fast_moving_variables)));
fast_macro_factors=zscore(table2array(fast_moving_variables))*coeff;
fast_macro_factors=fast_macro_factors(:,1:2);
fast_macro_factors=array2timetable(array2table(fast_macro_factors),'RowTimes', datetime(macro_data.Var1,'Format','yyyyMM'));
std(table2array(table2array(fast_macro_factors)));

plot(fast_macro_factors.Time,table2array(fast_macro_factors.Var1));
hold on 
plot(fast_macro_factors.Time,table2array(fast_macro_factors.Var2));

%let's reorder slow moving variables
slow_moving_variables=[slow_moving_variables(:,1) slow_moving_variables(:,48) slow_moving_variables(:,2:47) slow_moving_variables(:,49:65)];
%reordering fast moving variables
fast_moving_variables=[fast_moving_variables(:,1) fast_moving_variables(:,6) fast_moving_variables(:,2:5) fast_moving_variables(:,7:11)];

for i=1:size(slow_moving_variables,2)
pca_slow_factor_loadings(:,i)=inv(table2array(table2array(slow_macro_factors))'*table2array(table2array(slow_macro_factors)))*table2array(table2array(slow_macro_factors))'*table2array(slow_moving_variables(:,i));
end

for i=1:size(fast_moving_variables,2)
pca_fast_factor_loadings(:,i)=inv(table2array(table2array(fast_macro_factors))'*table2array(table2array(fast_macro_factors)))*table2array(table2array(fast_macro_factors))'*table2array(fast_moving_variables(:,i));
end

named_slow_macro_factors=table2array(table2array(slow_macro_factors))*pca_slow_factor_loadings(:,1:2);
corr(named_slow_macro_factors(:,1),named_slow_macro_factors(:,2));
named_slow_macro_factors=table2timetable([table(datetime(slow_moving_variables.Var1,'Format','yyyyMM')), array2table(named_slow_macro_factors)]);

named_fast_macro_factors=table2array(table2array(fast_macro_factors))*pca_fast_factor_loadings(:,1:2);
corr(named_fast_macro_factors(:,1),named_fast_macro_factors(:,2));
named_fast_macro_factors=table2timetable([table(datetime(slow_moving_variables.Var1,'Format','yyyyMM')), array2table(named_fast_macro_factors)]);


regression_factors=synchronize(named_slow_macro_factors,ffr,named_fast_macro_factors);
regression_factors=rmmissing(regression_factors);

full_sample=synchronize(named_slow_macro_factors,ffr,named_fast_macro_factors,slow_moving_variables,fast_moving_variables);
full_sample=rmmissing(full_sample);
named_slow_macro_factors=full_sample(:,1:2);
slow_moving_variables=full_sample(:,6:70);
constant=ones(size(named_slow_macro_factors,1),1);
for i=1:size(slow_moving_variables,2)
named_slow_factor_loadings(:,i)=inv([table2array(named_slow_macro_factors) constant]'*[table2array(named_slow_macro_factors) constant])*[table2array(named_slow_macro_factors) constant]'*[table2array(slow_moving_variables(:,i))];
end

for i=1:size(slow_moving_variables,2)
slow_variables_hat(:,i)=[table2array(named_slow_macro_factors) constant]*named_slow_factor_loadings(:,i);
end
slow_variables_residuals=table2array(slow_moving_variables)-slow_variables_hat;
slow_variables_residuals=array2timetable(slow_variables_residuals,'RowTimes',datetime(slow_moving_variables.Var1,'Format','yyyyMM'));

constant=ones(size(regression_factors,1),1);
fast_moving_variables=full_sample(:,71:81);
for i=1:size(fast_moving_variables,2)
named_fast_factor_loadings(:,i)=inv([table2array(regression_factors) constant]'*[table2array(regression_factors) constant])*[table2array(regression_factors) constant]'*table2array(fast_moving_variables(:,i));
end

for i=1:size(fast_moving_variables,2)
fast_variables_hat(:,i)=[table2array(regression_factors) constant]*named_fast_factor_loadings(:,i);
end
fast_variables_residuals=table2array(fast_moving_variables)-fast_variables_hat;
fast_variables_residuals=array2timetable(fast_variables_residuals,'RowTimes',datetime(fast_moving_variables.Var1,'Format','yyyyMM'));

factors=synchronize(regression_factors,betas);
factors=rmmissing(factors);
%Estimation of the state equation
Mdl = varm(8,1);
EstMdl = estimate(Mdl,factors);
var=summarize(EstMdl)
resid_covariance=var.Covariance;
H=chol(resid_covariance,'lower')*diag(diag(chol(resid_covariance,'lower')).^2)^(-1/2);
E = infer(EstMdl,factors);
response=irf(EstMdl);
armairf(EstMdl.AR,[],InnovCov=var.Covariance);
%Estimation of the IRFs
A_companion_form =[EstMdl.AR{1,1}]%; eye(8*11) zeros(8*11,8)];% 
eig_companion_form=eig(A_companion_form);
abs(eig_companion_form);

C(1:8,1:8)=eye(8);
years=5;
months=years*12;
for i=1:(months-1)
    bla=A_companion_form^i;
C(1:8,(8*i+1):(8*(i+1)))=bla(1:8,1:8);
end

inv_H=inv(H);
B(1:8,1:8)=inv_H;

for i=1:(months-1)
B(1:8,(8*i+1):(8*(i+1)))=C(1:8,(8*i+1):(8*(i+1)))*inv_H;
end

yields_loadings=[zeros(15,5), ns_factor_loadings];
full_named_slow_factor_loadings=[zeros(5,65); named_slow_factor_loadings];
full_named_fast_factor_loadings=[zeros(2,11); named_fast_factor_loadings];
complete_observation_equation_factor_loadings=[full_named_slow_factor_loadings full_named_fast_factor_loadings, yields_loadings'];
observed_variables_irf=complete_observation_equation_factor_loadings'*B;

plot(B(3,3:8:480));
plot(B(6,3:8:480));
plot(B(7,3:8:480));
plot(B(8,3:8:480));
plot(observed_variables_irf(79,3:8:480));
%Now we will do Yamamoto (2016)
factors_starting_values=factors(1,:);
state_equation_residuals=E(:,9);
state_equation_residuals=table2array(state_equation_residuals)-mean(table2array(state_equation_residuals));
n_bootstrap_samples=1000;
bootstrap_factors=zeros(size(factors,1),size(factors,2),n_bootstrap_samples);
for z=1:n_bootstrap_samples
bootstrap_sample_points=randsample((size(state_equation_residuals,1)),(size(state_equation_residuals,1)));
state_equation_bootstrap_residuals=state_equation_residuals(bootstrap_sample_points,:);
bootstrap_factors(1,:,z)=table2array(factors_starting_values);
estimated_AR=[EstMdl.AR{1,1}];
for i=1:(size(bootstrap_factors,1)-12)
    internal_sum=estimated_AR(1:8,1:8)*bootstrap_factors((i),:,z)';
   bootstrap_factors(i+12,:,z)= internal_sum+state_equation_bootstrap_residuals(i,:)';
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
bootstrap_observables(i,:,z)=[bootstrap_factors(i,1:2,z) 1 bootstrap_factors(i,3:8,z)]*complete_observation_equation_factor_loadings+observable_variables_bootstrap_residuals(i,:);
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
Mdl = varm(8,1);
EstMdl = estimate(Mdl,bootstrap_factors(:,:,z));
var=summarize(EstMdl);
resid_covariance=var.Covariance;
H=chol(resid_covariance,'lower')*diag(diag(chol(resid_covariance,'lower')).^2)^(-1/2);
E = infer(EstMdl,bootstrap_factors(:,:,z));
A_companion_form =[EstMdl.AR{1,1}];%; eye(8*11) zeros(8*11,8)];% 
eig_companion_form=eig(A_companion_form);
abs(eig_companion_form);

C(1:8,1:8)=eye(8);
years=5;
months=years*12;
for i=1:(months-1)
    bla=A_companion_form^i;
C(1:8,(8*i+1):(8*(i+1)))=bla(1:8,1:8);
end

inv_H=inv(H);
B_bootstrap(1:8,1:8,z)=inv_H;

for i=1:(months-1)
B_bootstrap(1:8,(8*i+1):(8*(i+1)),z)=C(1:8,(8*i+1):(8*(i+1)))*inv_H;
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

plot(B(3,3+(8*(0:59))));
hold on
plot(B(3,3+(8*(0:59)))+irf_standard_deviation*2);
hold on
plot(B(3,3+(8*(0:59)))-irf_standard_deviation*2);
title('Impulse Response Function: Interest Rate')
xlabel('Time period')
ylabel('Effect')

for i=0:59
irf_standard_deviation(i+1)=std(B_bootstrap(6,(3+8*i)+480*(1:999)));
end

plot(B(6,3+(8*(0:59))));
hold on
plot(B(6,3+(8*(0:59)))+irf_standard_deviation*2);
hold on
plot(B(6,3+(8*(0:59)))-irf_standard_deviation*2);
title('Impulse Response Function: Nelson-Siegel Yield Curve Level')
xlabel('Time period')
ylabel('Effect')

for i=0:59
irf_standard_deviation(i+1)=std(B_bootstrap(7,(3+8*i)+480*(1:999)));
end

plot(-1*B(7,3+(8*(0:59))));
hold on
plot(-1*B(7,3+(8*(0:59)))+irf_standard_deviation*2);
hold on
plot(-1*B(7,3+(8*(0:59)))-irf_standard_deviation*2);
title('Impulse Response Function: Nelson-Siegel Yield Curve Slope')
xlabel('Time period')
ylabel('Effect')

