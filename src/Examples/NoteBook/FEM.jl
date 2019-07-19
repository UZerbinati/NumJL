#Esercizi 6.24,6.25,6.26,7.2,7.4
include("ODEBasic.jl");
using PyPlot
q(x)=0*x+1;
f(x)=(1+pi^2)*sin(pi*x);
D=[0,1];
N=4;
y, x = FEM1H(q,f,D,N);
xx=D[1]:0.001:D[2];
figure()
title("Homogeneus FEM");
plot(xx,sin(pi*xx),label="Exact Solution");
plot(x,y,"*-",label="Numerical Solution");
legend(loc=2,borderaxespad=0)

q(x)=1/((1+x)^2);
f(x)=-1/((1+x)^3);
D=[0.0,1.0];
I=[1.0,0.5];
N=4;
figure()
y, x=FEM1(q,f,D,I,N);
xx=D[1]:0.001:D[2];
U =  [];
for p in xx
	U=push!(U,(1.0/(1.0+p)));
end
plot(xx,U,label="Exact Solution");
plot(x,y,"*-",label="Numerical Solution");
legend(loc=2,borderaxespad=0)
title("Dirichlet FEM");
show();
