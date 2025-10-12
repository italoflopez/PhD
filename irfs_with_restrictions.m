%IRFs with restrictions
fac_est_out_irf_restrictions = fac_est_out;
fac_est_out_irf_restrictions.varout.betahat([4,12,20], [6,7,8]) = 0;

irf_vdecomp_out_restrictions = dynamic_factor_irf_vdecomp(fac_est_out_irf_restrictions, est_par, decomp_par, bptcodevec);