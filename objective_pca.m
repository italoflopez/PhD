%Least Squares Objective Function
function value=objective_pca(x,y,z)
value=(size(x,2)*size(x,1))^(-1)*(sum((x-y*z)*(x-y*z')));
end