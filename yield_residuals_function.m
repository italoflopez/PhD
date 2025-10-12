function [estimated_yield_curve,yield_residuals]=yield_residuals_function(yield_data,ns_factor_loadings,betas)
estimated_yield_curve=ns_factor_loadings*table2array(betas)';
estimated_yield_curve=estimated_yield_curve';
yield_residuals=table2array(yield_data)-estimated_yield_curve;
yield_residuals=table2timetable([table(datetime(yield_data.Var1,'Format','yyyyMM')), array2table(yield_residuals)]);
estimated_yield_curve=array2timetable(array2table(estimated_yield_curve),'RowTimes', datetime(yield_data.Var1,'Format','yyyyMM'));
end