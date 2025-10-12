%Code for estimating the number of factors for slow moving variables
cd("C:\Users\Italo\OneDrive\Documents\PhD\PhD\First")
addpath('C:\Users\Italo\OneDrive\Documents\PhD\PhD\First\ddisk\matlab');
addpath('C:\Users\Italo\OneDrive\Documents\PhD\PhD\First\ddisk\m_utilities');
addpath('C:\Users\Italo\OneDrive\Documents\PhD\PhD');

%Now we will read the macro data
macro_data=readtable("stationary_data_for_macro_factors.csv");
macro_data=table2timetable([table(datetime(macro_data.Var1,'Format','yyyyMM')), macro_data(:,2:size(macro_data,2))]);

slow_moving_variables=macro_data(:,2:47);

est_par.smpl_par.calvec = 1:size(slow_moving_variables,1);
est_par.smpl_par.nper = 12;  % Monthly data
est_par.smpl_par.nfirst = [1,1]; % Index of the first observation
est_par.smpl_par.nlast = [length(est_par.smpl_par.calvec),length(est_par.smpl_par.calvec)];  % Last observation for estimation

est_par.fac_par.w = [];  % Matrix of observed factors (if applicable)
est_par.fac_par.lambda_constraints_est= 1;  % Constraint matrix
est_par.fac_par.nt_min = 50;  % Minimum number of observations for factor estimation
est_par.fac_par.tol = 1e-6;  % Precision of the estimate % Number of unobserved factors
est_par.fac_par.nfac.observed = 0;
est_par.fac_par.nfac.total = est_par.fac_par.nfac.unobserved + est_par.fac_par.nfac.observed;
est_par.n_uarlag = 4; %number of AR lags for uniquesess
est_par.lambda.nt_min = 100;

est_par.fac_par.nfac.unobserved = 4;

est_data = timetable2table(slow_moving_variables);
est_data = table2array(est_data(:, 2:end));



lsout = factor_estimation_ls(est_data, est_par);


bn_icp = bai_ng(lsout);

bn_icps=n_pca_factors(est_data,est_par,n_pca_max);

[min_num,min_idx] = min(bn_icps);