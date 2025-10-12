%Function for the number of pca factors


function [bn_icps]=n_pca_factors(est_data,est_par,n_pca_max)
bn_icps=zeros(n_pca_max,1);
for i = 1:n_pca_max
    est_par.fac_par.nfac.unobserved = i;
    est_par.fac_par.nfac.total = est_par.fac_par.nfac.unobserved + est_par.fac_par.nfac.observed;
    lsout = factor_estimation_ls(est_data, est_par);
    bn_icps(i,1)=bai_ng(lsout);
end
end