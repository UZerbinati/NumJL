using PyPlot
include("ODEBasic.jl");

N=40;
Eig=[];
NumEig=[];
q(x)= 0*x;
f(x)= 0*x;
D=[0,1];
I=[0,0];
for i in 1:N
	Eig=push!(Eig,((pi*i)/1)^2);
	
end
# -u''(x)+q(x)u(x)=f(x);
y,x,h, A = FDirichlet(q,f,D,I,N);
NumEig, NumEigV = eigs(A; nev=N);

NumEig = sort(NumEig);

println(Eig);
println(NumEig);

figure();
title("EigenValues");
plot(Eig,"bo",label="Exact Egenvalues");
plot(NumEig,"r*",label="Eienvalue of Finite Difference Discretization");
legend(loc=2,borderaxespad=0)

K = [2,5,20,39];
t=0:0.001:1;
for l in K
	println(l);
	figure();
	u(x)= sin(((pi*l)/1)*x);
	y,x,h, A = FDirichlet(q,f,D,I,N);
	NumEig, NumEigV = eigs(A; nev=N);
	println("TEST");
	println(u(t));
	println(size(A));
	println(length(0:(1/(N-1)):1));
	plot(t,u(t));
	v=NumEigV[:,40-l];
	norma = norm(u(t),2)*sign(u(0));
	print("norma",norma);
	norma = norma/(norm(v,2)*sign(v[1]));
	plot(0:(1/(N-1)):1,v*norma,"-o");
end

show();

