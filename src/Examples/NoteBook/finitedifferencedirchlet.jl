include("ODEBasic.jl");
using PyPlot;
u(x)= cosh(x)/cosh(1);
q(x)=0*x+1;
f(x)=0*x;
y,x =FDirichlet(q,f,[-1.0,1.0],[1.0,1.0],20);
figure();
plot(x,y,"b*",label="Finite Difference Approximation");
plot(x,u(x),"r-",label="Exact Solution");
legend(loc=2,borderaxespad=0);
show();
figure();
H=[];
U2=[];
U1=[];
Uinf=[];
C = [];
for N in 2:2:18
	y,x,h, c =FDirichlet(q,f,[-1.0,1.0],[1.0,1.0],2^N);
	H=push!(H,h);
	U1=push!(U1,norm(y-u(x),1));
	U2=push!(U2,norm(y-u(x),2));
	Uinf=push!(Uinf,norm(y-u(x),Inf));
end
println("End Computation");
loglog(H,U1,label="Approximation Error Norm 1");
loglog(H,U2,label="Approximation Error Norm 2");
loglog(H,Uinf,label="Approximation Error Norm Inf");
legend(loc=2,borderaxespad=0)
show();
