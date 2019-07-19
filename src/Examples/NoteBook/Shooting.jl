using PyPlot;
include("ODEBasic.jl");
du1(x,y,t)=y;
du2(x,y,t)=9.8;
y=shootingBVP(du1,du2,[0.0,5.0],[0.0,40.0],[1.0,50.0],0.1);
figure()
plot(t,y[1]);
figure()
plot(t,y[2]);
show()
