%lag order aic
function [BICs]=lag_order_selection_BIC(data,pmax)
BICs=zeros(pmax,1);
for i = 1:pmax
Mdl = varm(size(data,2),i);
EstMdl = estimate(Mdl,data(1:(size(data,1)-pmax+i-1),:));%
var=summarize(EstMdl);
BICs(i,1)=var.BIC;
end