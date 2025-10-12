function [betas,ns_factor_loadings] = betas_function(yield_data)

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
ns_factor_loadings;
end


