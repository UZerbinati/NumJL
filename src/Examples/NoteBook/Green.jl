include("ODEBasic.jl");
using PyPlot;
u(x)= cosh(x)/cosh(1);
q(x)=0*x+1;
f(x)=0*x;
y,x,h,A =FD(q,f,[-1.0,1.0],[1.0,1.0],30);
surf(inv(full(A)),cmap=ColorMap("plasma"));
show();
