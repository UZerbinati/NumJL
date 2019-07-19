#Example: "Lotka-Volterra"
using PyPlot;
include("ODEBasic.jl");
du1(x,y,t)=8*x-2*x*y;
du2(x,y,t)=-5*y+3*x*y;
DU(X,t)=[du1(X[1],X[2],t),du2(X[3],X[4],t)];
y,t=eulerND(DU,[0.0,20.0],0.0001,[3.0,4.0]);
figure();
plot(t,y[2]);
plot(t,y[1]);
show();
