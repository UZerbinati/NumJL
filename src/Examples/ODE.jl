function LotkaVolterraExample()
	du1(x,y,t)=8*x-2*x*y;
	du2(x,y,t)=-5*y+3*x*y;
	y, t=euler2D(du1,du2,[0.0,20.0],0.000001,[3.0,4.0])
	plot(t,y[2]);
	plot(t,y[1]);
end
function VanderPolExample()
	μ=1;
	du1(x,y,t)=μ*(x-(1/3)*x^3-y);
	du2(x,y,t)=(1/μ)*x;
	y, t=euler2D(du1,du2,[0.0,100.0],0.001,[3.0,4.0])
	figure();
	plot(t,y[2]);
	plot(t,y[1]);
	figure();
	plot(y[1],y[2]);
end
function StabilityRegionEuler()
	N = 10000;
	D = [-5.0,5.0];
	Z = ones(ComplexF64,N,N);
	K = zeros(N,N);
	h=(D[2]-D[1])/(N+1);
	Y,X = ndgrid(D[1]:h:D[2],D[1]:h:D[2]);
	for i in 1:N
		for j in 1:N
			while abs(Z[i,j])> 0.0000001 && K[i,j] < 2
				du(x,t)=(X[i,j]+Y[i,j]*1im)*x;
				y, t = eulerMethod(du,[0.0,50.0],0.9*(K[i,j]+1),1.0+0.0im);
				Z[i,j]=y[end];
				#semilogy(t,abs.(y));
				#readline()
				K[i,j]=K[i,j]+1;
			end
		end
		println(i/N);
	end
	imshow(K, cmap="hot", interpolation="nearest",extent=[-5.0,5.0,-5.0,5-0]);
	colorbar();
	title("Euler Absolute Stability Region")
end
function StabilityRegionImplicitEuler()
	N = 1000;
	D = [-5.0,5.0];
	Z = zeros(ComplexF64,N,N);
	K = zeros(N,N);
	h=(D[2]-D[1])/(N+1);
	Y,X = ndgrid(D[1]:h:D[2],D[1]:h:D[2]);
	for i in 1:N
		for j in 1:N
			while abs(Z[i,j])< 0.000001 && K[i,j] < 2
				du(x,t)=(X[i,j]+Y[i,j]*1im)*x;
				ddu(x,t)=(X[i,j]+Y[i,j]*1im);
				y, t = complexEulerImplicitMethod(du,ddu,[0.0,250.0],0.9*(K[i,j]+1),1.0+0.0im);
				Z[i,j]=y[end];
				K[i,j]=K[i,j]+1;
			end
		end
		println(i/N);
	end
	imshow(K, cmap="hot", interpolation="nearest",extent=[-5.0,5.0,-5.0,5-0]);
	colorbar();
	title("Implicit Euler Absolute Stability Region");
end
