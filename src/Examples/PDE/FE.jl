function FExample()
	q(x)=1/((1+x)^2);
	f(x)=-1/((1+x)^3);
	D=[0.0,1.0];
	Ini=[1.0,0.5];
	N=4;
	figure()
	y, x=FEM1(q,f,D,Ini,N);
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
end
