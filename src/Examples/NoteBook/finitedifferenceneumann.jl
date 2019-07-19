include("ODEBasic.jl");
using PyPlot
u(x)= 1/x;
q(x)=1/(x^2);
f(x)=(-1)*x^(-3);
y1,x1 =FDNeumannSlow(q,f,[1.0,2.0],[-1.0,-1/4],40);
y2,x2 =FDNeumann(q,f,[1.0,2.0],[-1.0,-1/4],40);
figure();
plot(x1,y1,"b*",label="Finite Difference Approximation Linear Method");
plot(x2,y2,"g*",label="Finite Difference Approximation Quadratic Method");
ex=[];
for i in 1:length(x1);
	ex=push!(ex,u(x1[i]));
end
plot(x1,ex,"r-",label="Exact Solution");
legend(loc=2,borderaxespad=0);
show();
