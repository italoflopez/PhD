function slow_macro_factors = factors_function(variables,number_of_factors)

[coeff,score,latent,tsquared,explained,mu] = pca(zscore(table2array(variables)));
slow_macro_factors=zscore(table2array(variables))*coeff;
number_of_factors=number_of_factors;
slow_macro_factors=slow_macro_factors(:,1:number_of_factors);
slow_macro_factors=array2timetable(array2table(slow_macro_factors),'RowTimes', datetime(variables.Var1,'Format','yyyyMM'));
slow_macro_factors
end