%lag order aic
function [AICs]=lag_order_selection(data,pmax)
AICs=zeros(pmax,1);
for i = 1:pmax
Mdl = varm(size(data,2),i);
EstMdl = estimate(Mdl,data(1:(size(data,1)-pmax+i-1),:));%
var=summarize(EstMdl)
AICs(i,1)=var.AIC;
end
