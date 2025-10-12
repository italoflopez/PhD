%function to calculate true number of factors
function [n_factors]=bai_ng_n_factors(data,max_n_factors)
zdata=zscore(table2array(data));
IC_r_vector=zeros(max_n_factors,1);
[coeff,score,latent,tsquared,explained,mu] = pca(zdata);
for i=1:max_n_factors
factors=zdata*coeff;
factors=factors(:,1:i)
matrix_inverted=inv(factors'*factors);
clear pca_factor_loadings
for j=1:size(zdata,2)
pca_factor_loadings(:,j)=matrix_inverted*factors'*zdata(:,j);
end
clear residuals
for k=1:size(zdata,2)
residuals(:,k)=zdata(:,k)-factors*pca_factor_loadings(:,k);
end
V_r=sum(sum(residuals'*residuals))/(size(residuals,1)*size(residuals,2))
IC_2=min(size(zdata,2),size(zdata,1))*((size(zdata,2)+size(zdata,1))/(size(zdata,2)*size(zdata,1)))*log(min(size(zdata,2),size(zdata,1)))
IC_r=V_r+i*IC_2
IC_r_vector(i)=IC_r
end
[min_num,min_idx] = min(IC_r_vector);
n_factors=min_idx