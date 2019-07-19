function finiteDifferenceExample()
	#Defining the function we are going to apply finite difference to.
	f(t)=cos(t)-sin(exp(t));
	#Defining the point where we are going to evaluate the derivative
	x=1;
	banchmark=-sin(1)-exp(1)*cos(exp(1));
	println("Exact value in x, for the derivative:", banchmark);
	H=[];
	DF=[];#Vector of farward finite difference
	DFC=[];#Vector of central finite difference
	DF4=[];#Vector of 4th order quadrature finite difference
	ComplexStepDiff=[];#Vector for complex step differentiation
	#Defining vector cotaing the difference between the computed derivative by each method and the banchmark.
	ErrorDF=[];
	ErrorDFC=[];
	ErrorDF4=[];
	ErrorComplexStepDiff=[];
	h=0.5;
	#Computing each finite difference and each error, populating the vector defined above.
	while (h>eps())
		println("Length of finite difference: ",h)
		H=append!(H,h);
		y=(f(x+h)-f(x))/h;
		yc=(f(x+h)-f(x-h))/(2*h);
		y4=(1/h)*((1/12)*f(x-2*h)-(2/3)*f(x-h)+(2/3)*f(x+h)-(1/12)*f(x+2*h));
		yi=imag(f(x+h*1im)/h);
		DF=append!(DF,y);
		DFC=append!(DFC,yc);
		DF4=append!(DF4,y4);
		ErrorDF=append!(ErrorDF,abs(banchmark-y));
		ErrorDFC=append!(ErrorDFC,abs(banchmark-yc));
		ErrorDF4=append!(ErrorDF4,abs(banchmark-y4));
		ErrorComplexStepDiff=append!(ErrorComplexStepDiff,abs(banchmark-yi));
		h=h*0.5;
	end
	println("Vector of finite differences length",H);
	figure()
	title("Finite Difference");
	loglog(H,ErrorDF,label="Farward Finite Difference");
	loglog(H,ErrorDFC,label="Central Finite Difference");
	loglog(H,ErrorDF4,label="4th Order Quadrature Finite Difference");
	loglog(H,ErrorComplexStepDiff,"r--",label="Complex Step Derivative");
	grid("on")
	legend(loc=2,borderaxespad=0)
	show();
end
function FDirichletExample()
    u(x)= cosh(x)/cosh(1);
    q(x)=0*x+1;
    f(x)=0*x;
    y,x =FDirichlet(q,f,[-1.0,1.0],[1.0,1.0],20);
	#println(cos.(collect(x)))
	figure();
    plot(x,y,"b*",label="Finite Difference Approximation");
    plot(x,u.(x),"r-",label="Exact Solution");
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
    	U1=push!(U1,norm(y-u.(x),1));
    	U2=push!(U2,norm(y-u.(x),2));
    	Uinf=push!(Uinf,norm(y-u.(x),Inf));
    end
    println("End Computation");
	#println(size(x))
	#println(size(y))
    loglog(H,U1,label="Approximation Error Norm 1");
    loglog(H,U2,label="Approximation Error Norm 2");
    loglog(H,Uinf,label="Approximation Error Norm Inf");
    legend(loc=2,borderaxespad=0)
    show();
end
function FDNeumannExample()
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
end
function FDScalarProductExample()
	D = [0,1];
	f(x)=(x-1)*x;
	df(x)=2*x-1;
	Error = [];
	for N in 10:10:30
		h = (D[2]-D[1])/(N+1);
		x = D[1]:h:D[2];
		U = f.(x);
		println(U);
		d1=ones(1,N-1);
		d1=convert(Array{Float64}, d1);
		d2=-2*ones(1,N);
		d2=convert(Array{Float64}, d2);
		A=spdiagm(-1=>(1)*d1[1,:],0=>d2[1,:],1=>(1)*d1[1,:]);
		A = (1/(h^2))*A;
		println(Matrix(A))
		println(size(A))
		y = A\U[2:end-1];
		Uh = [0.0];
		Uh = append!(Uh,y);
		Uh = push!(Uh,0.0);
		println(Uh)
		Error = push!(Error,norm(df.(x)-Uh,2));
		figure()
		plot(U)
		plot(df.(x))
		plot(Uh,"r*");
		readline();
	end
	figure()
	plot(10:10:30,Error);
end
