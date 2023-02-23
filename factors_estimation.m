%factors estimation
function [factors,factor_loadings]=factors_estimation(data,number_of_factors,restrictions)
x=data;
n=size(data,2);
r=number_of_factors;
objective_pca=@(y,z) trace(size(x,2)^(-1)*(sum((x-y*z)*(x-y*z)')));
prob.Objective=trace(size(x,2)^(-1)*(sum((x-y*z)*(x-y*z)')));
w0.y=1;
w0.z=[1 1]
