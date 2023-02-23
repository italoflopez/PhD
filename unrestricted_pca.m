%unrestricted pca
function value=unrestricted_pca(X)
r=1;
N=2;
T=100;
x_t=normrnd(0,1,N,T);
X=sym('X',[N*r+T*r,1]);factor_loadings=X(1:(N*r));
factors=X((N*r+1):end);
value=(size(x_t,2)*size(x_t,1))^(-1)*(sum((x_t-factor_loadings*factors').^2,'all'));
end
