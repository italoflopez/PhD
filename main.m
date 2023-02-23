objective=@(x) x(1)*x(4)*(x(1)+x(2)+x(3))+x(3);
x0=[1,5,5,1];
lb=ones(4);
ub=5*ones(4);
A=[];
B=[];
Aeq=[];
Beq=[];
nonlincon=@nlcon;

x=fmincon(objective,x0,A,B,Aeq,Beq,lb,ub,nonlincon);