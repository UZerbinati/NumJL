include("ODEBasic.jl");
using PyPlot
f(x)=sin(x);
df(x)=cos(x);
N=100;
I=[0.0,2*pi];
D=[pi/3,pi/3];
h=(I[2]-I[1])/(N+1);
println("Passo: ",h);
x=I[1]:h:I[2];
u=[]
for i in 2:N+1
	u=push!(u,f(x[i]/2)+pi/3)
end
println(u');
#u=zeros(1,N);
#println(u);
x, y = FDNonLinearDirchletAutonomus(f,df,I,D,N,u');
println(x);
println(length(x));
println(y);
println(length(y));
plot(x,y)
show();
