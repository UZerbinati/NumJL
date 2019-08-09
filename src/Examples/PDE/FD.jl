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
function Interpoler(z,w,X,Y,Z)
	m = findall(x -> x == z,X);
	n = findall(x -> x == w,Y);
	#println(Z[m,n])
	return norm(Z[m,n],Inf)
end
function spacer(X,N)
	Y=[];
	for i in 1:size(X)[1]-1
		Y = append!(Y,range(X[i],stop=X[i+1],length=N));
	end
	return Y;
end
function oddDomainExample(D)
	PN = [];
	PC = [];
	
	N = 60;
	r=10; #Stepsize for solution
	Delta = 0.01*(1/sqrt(2)); #Time step	
	X = range(-3,stop=3,length=N);
	Y = range(3,stop=-3,length=N);	
	X2 = spacer(X,r);
	Y2 = spacer(Y,r); 
	Ini = zeros(size(X2)[1],size(Y2)[1]);
	Ini2 = zeros(size(X2)[1],size(Y2)[1]);
	Int= zeros(size(X2)[1],size(Y2)[1]);

	
	μx=1.5; μy=-1.5; σ=0.1;
        f(x,y)=3*exp((-(x-μx)^2)/(2*σ^2)) * exp((-(y-μy)^2)/(2*σ^2));
	df(x,y) = 0*x*y;
	Gamma(x,y) = 0*x*y;
	
	# [INFO] Showing the domain
	Z = ones(N,N)
	for i in 1:size(X)[1]
		for j in 1:size(Y)[1]
			Z[i,j]=D(X[i],Y[j]);
		end
	end
	show(Z);
	matshow(Z,origin="upper",extent=[-3,3,-3,3]);
	#colorbar();
	#figure()
	title("Meshing...");

	# [INFO] Some sort of Meshing

	for i in 1:size(X)[1]
		for j in 1:size(Y)[1]
			if Z[i,j] == 1
				if sum(Z[1:i,j]) == 1
					plot(X[i],Y[j],"g*");
					Z[i,j]=2;
					PN = push!(PN,[i,j]);
					PC = push!(PC,[X[i],Y[j]]);
                                elseif sum(Z[i,j:size(Z)[2]])==1
                                        plot(X[i],Y[j],"g*");
                                        Z[i,j]=2;
					PN = push!(PN,[i,j]);
                                        PC = push!(PC,[X[i],Y[j]]);
				elseif sum(Z[i,1:j]) == 1
                                        plot(X[i],Y[j],"g*");
                                        Z[i,j]=2;
					PN = push!(PN,[i,j]);
                                        PC = push!(PC,[X[i],Y[j]]);
				elseif sum(Z[i:size(Z)[1],j])==1
					plot(X[i],Y[j],"g*");
					Z[i,j]=2;
					PN = push!(PN,[i,j]);
					PC = push!(PC,[X[i],Y[j]]);
				else
					plot(X[i],Y[j],"r*");
				end
			end
		end
	end

	show(Z);

	# [INFO] Builing the square where we are solving the equation
	I = 0;
	println("--------------------!WHERE ARE WE SOLWING!--------------------");
	while true
		S = zeros(size(X2)[1],size(Y2)[1]);
		println("----------|Solving Step|-----------");	
		I = I+1;
		for i in 1:N-1
			for j in 1:N-1
				if Z[i,j] == 1
					if Z[i+1,j] == 1
						if Z[i,j+1] == 1
							if Z[i+1,j+1]==1
							#   *....*
							#   :    :
							#   *....*
								
								println("Solving only on inside, red O",[X[i],X[i+1]],[Y[j],Y[j+1]]);
								#plot(X2[(r*i)-r+1],Y2[(r*j)-r+1],"ms");
								H = X[2]-X[1];
								h = X2[2]-X2[1];
								#plot(X2[r*i-r]+(H/2),Y2[r*j-r]-(H/2),"ro");
								for m in 1:r
									for n in 1:r
										k1 = (r*i-r)+m;
										k2 = (r*j-r)+n;
										if I == 1
											Ini[k1,k2] = f(X2[k1],Y2[k2]);	
											Ini2[k1,k2] = df(X2[k1],Y2[k2]);
										elseif I == 2
											S[k1,k2]=f(X2[k1],Y2[k2])+Delta*df(X2[k1],Y2[k2])+0.5*(Delta^2/(h^2))*(Ini[k1+1,k2]-2*Ini[k1,k2]+Ini[k1-1,k2]+Ini[k1,k2-1]-2*Ini[k1,k2]+Ini[k1,k2+1]);
										elseif I>2
											#println("SBAM");
											#println(size((Ini[end-1])));
											S[k1,k2]=(Delta^2/(h^2))*(Ini[end][k1+1,k2]-2*Ini[end][k1,k2]+Ini[end][k1-1,k2]+Ini[end][k1,k2-1]-2*Ini[end][k1,k2]+Ini[end][k1,k2+1])+2*Ini[end][k1,k2]-Ini[end-1][k1,k2];
										end
									end
								end
							
					
							elseif Z[i+1,j+1]==2
							#   *....*
							#   :    :
							#   *....o
								println("Solving only 1 point, white O");
								H = X[2]-X[1];
								h = X2[2]-X2[1];
        							plot(X2[r*i-r]+(H/2),Y2[r*j-r]-(H/2),"wo");
								for m in 1:r
									for n in 1:r
										if I == 1
											k1 = (r*i)-r+m;
											k2 = (r*j)-r+n;
											Ini[k1,k2] = f(X2[k1],Y2[k2]);	
											Ini2[k1,k2] = df(X2[k1],Y2[k2]);
										elseif I == 2
											S[k1,k2]=f(X2[k1],Y2[k2])+Delta*df(X2[k1],Y2[k2])+(Delta^2/(h^2))*(Ini[k1+1,k2]-2*Ini[k1,k2]+Ini[k1-1,k2]+Ini[k1,k2-1]-2*Ini[k1,k2]+Ini[k1,k2+1]);

										elseif I>2
											#println("SBAM");
											#println(size((Ini[end-1])));
											S[k1,k2]=(Delta^2/(h^2))*(Ini[end][k1+1,k2]-2*Ini[end][k1,k2]+Ini[end][k1-1,k2]+Ini[end][k1,k2-1]-2*Ini[end][k1,k2]+Ini[end][k1,k2+1])+2*Ini[end][k1,k2]-Ini[end-1][k1,k2];
										end
									end
								end
							end
						elseif Z[i,j+1] == 2
							if Z[i+1,j+1]==1
							#   *....*
							#   :    :
							#   o....*
								println("Solving with 1 points, white O");
								H = X[2]-X[1];
                                                                h = X2[2]-X2[1];
                                                                plot(X2[r*i-r]+(H/2),Y2[r*j-r]-(H/2),"wo");
								for m in 1:r
									for n in 1:r
										k1 = (r*i)-r+m;
										k2 = (r*j)-r+n; if I == 1
											Ini[k1,k2] = f(X2[k1],Y2[k2]);	
											Ini2[k1,k2] = df(X2[k1],Y2[k2]);
										elseif I == 2
											S[k1,k2]=f(X2[k1],Y2[k2])+Delta*df(X2[k1],Y2[k2])+(Delta^2/(h^2))*(Ini[k1+1,k2]-2*Ini[k1,k2]+Ini[k1-1,k2]+Ini[k1,k2-1]-2*Ini[k1,k2]+Ini[k1,k2+1]);
										elseif I>2
											#println("SBAM");
											#println(size((Ini[end-1])));
											S[k1,k2]=(Delta^2/(h^2))*(Ini[end][k1+1,k2]-2*Ini[end][k1,k2]+Ini[end][k1-1,k2]+Ini[end][k1,k2-1]-2*Ini[end][k1,k2]+Ini[end][k1,k2-1])+2*Ini[end][k1,k2]-Ini[end-1][k1,k2];
										end
									end
								end
							elseif Z[i+1,j+1] == 2
							#   *....*
							#   :    :
							#   o____o
								println("Solving with 1 border, green O");
								H = X[2]-X[1];
                                                                h = X2[2]-X2[1];
                                                                plot(X2[r*i-r]+(H/2),Y2[r*j-r]-(H/2),"go");
								for m in 1:r
									for n in 1:r
										k1 = (r*i)-r+m;
										k2 = (r*j)-r+n;
										if I == 1
											Ini[k1,k2] = f(X2[k1],Y2[k2]);	
											#Ini2[k1,k2] = df(X2[k1],Y2[k2]);
										elseif I == 2
											S[k1,k2]=f(X2[k1],Y2[k2])+Delta*df(X2[k1],Y2[k2])+(Delta^2/(h^2))*(Ini[k1+1,k2]-2*Ini[k1,k2]+Ini[k1-1,k2]+Ini[k1,k2-1]-2*Ini[k1,k2]+Ini[k1,k2+1]);
										elseif I>2
											#println("SBAM");
											#println(size((Ini[end-1]));
											S[k1,k2]=(Delta^2/(h^2))*(Ini[end][k1+1,k2]-2*Ini[end][k1,k2]+Ini[end][k1-1,k2]+Ini[end][k1,k2-1]-2*Ini[end][k1,k2]+Ini[end][k1,k2+1])+2*Ini[end][k1,k2]-Ini[end-1][k1,k2];
										end
									end
								end
								if I > 1 
									YSLength = r;
									for m in 1:r
										k1 = (r*i)-r+m;
                                                                                S[k1,(r*j)-r+YSLength]=Gamma(X2[k1],Y2[(r*j)-r+YSLength])
										#plot(X2[(r*i)-r+m],Y2[(r*j-r)+YSLength],"y*");
									end
								end
							end
						end
					elseif Z[i+1,j]==2
						#println("OK",Z[i,j+1]);
						if Z[i,j+1] == 1
							if Z[i+1,j+1]==1
							#   *....o
							#   :    :
							#   *....*
								println("Solving with 1 point, white O");
								H = X[2]-X[1];
								h = X2[2]-X2[1];
                                                                plot(X2[r*i-r]+(H/2),Y2[r*j-r]-(H/2),"wo");
								for m in 1:r
									for n in 1:r
										k1 = (r*i)-r+m;	
	                                                                        k2 = (r*j)-r+n;

										if I == 1
											Ini[k1,k2] = f(X2[k1],Y2[k2]);	
											Ini2[k1,k2] = df(X2[k1],Y2[k2]);
										elseif I == 2
											S[k1,k2]=f(X2[k1],Y2[k2])+Delta*df(X2[k1],Y2[k2])+(Delta^2/(h^2))*(Ini[k1+1,k2]-2*Ini[k1,k2]+Ini[k1-1,k2]+Ini[k1,k2-1]-2*Ini[k1,k2]+Ini[k1,k2+1]);
										elseif I>2
											#println("SBAM");
											#println(size((Ini[end-1])));
											S[k1,k2]=(Delta^2/(h^2))*(Ini[end][k1+1,k2]-2*Ini[end][k1,k2]+Ini[end][k1-1,k2]+Ini[end][k1,k2-1]-2*Ini[end][k1,k2]+Ini[end][k1,k2+1])+2*Ini[end][k1,k2]-Ini[end-1][k1,k2];
										end
									end
								end
							elseif Z[i+1,j+1] == 2
							#   *....o
							#   :    |
							#   *....o
								println("Solving with 1 border");
								H = X[2]-X[1];
                                                                h = X2[2]-X2[1];
                                                                plot(X2[r*i-r]+(H/2),Y2[r*j-r]-(H/2),"go");
								for m in 1:r
									for n in 1:r
										k1 = (r*i)-r+m;
										k2 = (r*j)-r+n;
										if I == 1
											Ini[k1,k2] = f(X2[k1],Y2[k2]);	
											Ini2[k1,k2] = df(X2[k1],Y2[k2]);
										elseif I == 2
											S[k1,k2]=f(X2[k1],Y2[k2])+Delta*df(X2[k1],Y2[k2])+(Delta^2/(h^2))*(Ini[k1+1,k2]-2*Ini[k1,k2]+Ini[k1-1,k2]+Ini[k1,k2-1]-2*Ini[k1,k2]+Ini[k1,k2+1]);
										elseif I>2
											#println("SBAM");
											#println(size((Ini[end-1])));
											S[k1,k2]=(Delta^2/(h^2))*(Ini[end][k1+1,k2]-2*Ini[end][k1,k2]+Ini[end][k1-1,k2]+Ini[end][k1,k2-1]-2*Ini[end][k1,k2]+Ini[end][k1,k2+1])+2*Ini[end][k1,k2]-Ini[end-1][k1,k2];
										end
									end
								end
								if I > 1 
									XSLength = r;
									for m in 1:r
                                                                        	k2 = (r*j)-r+m;
										S[(r*i-r)+XSLength,k2]=Gamma(X2[(r*j)-r+XSLength],Y2[k2])
									end
								end
							end
						elseif Z[i,j+1] == 2
							if Z[i+1,j+1] == 2
							#   *....o
							#   :    |
							#   o____o
								println("Solving with 2 border, blue O");
								H = X[2]-X[1];
                                                                h = X2[2]-X2[1];
                                                                plot(X2[r*i-r]+(H/2),Y2[r*j-r]-(H/2),"go");
								for m in 1:r
									for n in 1:r
									k1 = (r*i)-r+m;
									k2 = (r*j)-r+n; 	
									if I == 1
											Ini[k1,k2] = f(X2[k1],Y2[k2]);	
											Ini2[k1,k2] = df(X2[k1],Y2[k2]);
										elseif I == 2
											S[k1,k2]=f(X2[k1],Y2[k2])+Delta*df(X2[k1],Y2[k2])+(Delta^2/(h^2))*(Ini[k1+1,k2]-2*Ini[k1,k2]+Ini[k1-1,k2]+Ini[k1,k2-1]-2*Ini[k1,k2]+Ini[k1,k2+1]);
										elseif I>2
											#println("SBAM");
											#println(size((Ini[end-1])));
											S[k1,k2]=(Delta^2/(h^2))*(Ini[end][k1+1,k2]-2*Ini[end][k1,k2]+Ini[end][k1-1,k2]+Ini[end][k1,k2-1]-2*Ini[end][k1,k2]+Ini[end][k1,k2+1])+2*Ini[end][k1,k2]-Ini[end-1][k1,k2];
										end
									end
								end
								if I > 1 
									XSLength = r;
									YSLength = r;
									for m in 1:r
										k2 = r*j-r+m;
										S[(r*i)-r+XSLength,k2]=Gamma(X2[(r*j)-r+XSLength],Y2[(r*j)-r+YSLength])
									end
									for m in 1:r
										k1=r*i-r+m;
										S[k1,(r*j)-r+YSLength]=Gamma(X2[(r*j)-r+XSLength],Y2[(r*j)-r+YSLength]) 	
									end
									
								end
							end
							if Z[i+1,j+1] == 1
							#   *....o
							#   :    :
							#   o----*
								println("Solving with 2 points, white O");
								H = X[2]-X[1];
                                                                h = X2[2]-X2[1];
                                                                plot(X2[r*i-r]+(H/2),Y2[r*j-r]-(H/2),"wo");
								for m in 1:r
									for n in 1:r
										k1 = (r*i)-r+m;
										k2 = (r*j)-r+n;
										if I == 1
											Ini[k1,k2] = f(X2[k1],Y2[k2]);	
											Ini2[k1,k2] = df(X2[k1],Y2[k2]);
										elseif I == 2
											S[k1,k2]=f(X2[k1],Y2[k2])+Delta*df(X2[k1],Y2[k2])+(Delta^2/(h^2))*(Ini[k1+1,k2]-2*Ini[k1,k2]+Ini[k1-1,k2]+Ini[k1,k2-1]-2*Ini[k1,k2]+Ini[k1,k2+1]);
										elseif I>2
											#println("SBAM");
											#println(size((Ini[end-1])));
											S[k1,k2]=(Delta^2/(h^2))*(Ini[end][k1+1,k2]-2*Ini[end][k1,k2]+Ini[end][k1-1,k2]+Ini[end][k1,k2-1]-2*Ini[end][k1,k2]+Ini[end][k1,k2+1])+2*Ini[end][k1,k2]-Ini[end-1][k1,k2];
										end
									end
								end
							end
						end

					end
				 elseif Z[i,j] == 2
					if Z[i+1,j] == 1
						if Z[i,j+1] == 1
							if Z[i+1,j+1]==1
							#   o....*
							#   :    :
							#   *....*
								println("Solving with 1 points white O");
								H = X[2]-X[1];
                                                                h = X2[2]-X2[1];
                                                                plot(X2[r*i-r]+(H/2),Y2[r*j-r]-(H/2),"wo");
								for m in 1:r
									for n in 1:r
										k1 = (r*i)-r+m;
										k2 = (r*j)-r+n;
										if I == 1
											Ini[k1,k2] = f(X2[k1],Y2[k2]);	
											Ini2[k1,k2] = df(X2[k1],Y2[k2]);
										elseif I == 2
											S[k1,k2]=f(X2[k1],Y2[k2])+Delta*df(X2[k1],Y2[k2])+(Delta^2/(h^2))*(Ini[k1+1,k2]-2*Ini[k1,k2]+Ini[k1-1,k2]+Ini[k1,k2-1]-2*Ini[k1,k2]+Ini[k1,k2+1]);
										elseif I>2
											#println("SBAM");
											#println(size((Ini[end-1])));
											S[k1,k2]=(Delta^2/(h^2))*(Ini[end][k1+1,k2]-2*Ini[end][k1,k2]+Ini[end][k1-1,k2]+Ini[end][k1,k2-1]-2*Ini[end][k1,k2]+Ini[end][k1,k2+1])+2*Ini[end][k1,k2]-Ini[end-1][k1,k2];
										end
									end
								end
							elseif Z[i+1,j+1] == 2
							#   o....*
							#   :    :
							#   *....o
								println("Solving with 2 points, white O");
								H = X[2]-X[1];
                                                                h = X2[2]-X2[1];
                                                                plot(X2[r*i-r]+(H/2),Y2[r*j-r]-(H/2),"wo");
								for m in 1:r
									for n in 1:r
										k1 = (r*i)-r+m;
										k2 = (r*j)-r+n; 
										if I == 1
											Ini[k1,k2] = f(X2[k1],Y2[k2]);	
											Ini2[k1,k2] = df(X2[k1],Y2[k2]);
										elseif I == 2
											S[k1,k2]=f(X2[k1],Y2[k2])+Delta*df(X2[k1],Y2[k2])+(Delta^2/(h^2))*(Ini[k1+1,k2]-2*Ini[k1,k2]+Ini[k1-1,k2]+Ini[k1,k2-1]-2*Ini[k1,k2]+Ini[k1,k2+1]);
										elseif I>2
											#println("SBAM");
											#println(size((Ini[end-1])));
											S[k1,k2]=(Delta^2/(h^2))*(Ini[end][k1+1,k2]-2*Ini[end][k1,k2]+Ini[end][k1-1,k2]+Ini[end][k1,k2-1]-2*Ini[end][k1,k2]+Ini[end][k1,k2+1])+2*Ini[end][k1,k2]-Ini[end-1][k1,k2];
										end
									end
								end
							end
						elseif Z[i,j+1] == 2
							if Z[i+1,j+1] == 2
							#   o....*
							#   |    :
							#   o____o
								println("Solving with 2 border, blue O");
								H = X[2]-X[1];
                                                                h = X2[2]-X2[1];
                                                                plot(X2[r*i-r]+(H/2),Y2[r*j-r]-(H/2),"go");
								for m in 1:r
									for n in 1:r
										k1 = (r*i)-r+m;
										k2 = (r*j)-r+n; 
										if I == 1
											Ini[k1,k2] = f(X2[k1],Y2[k2]);	
											Ini2[k1,k2] = df(X2[k1],Y2[k2]);
										elseif I == 2
											S[k1,k2]=f(X2[k1],Y2[k2])+Delta*df(X2[k1],Y2[k2])+(Delta^2/(h^2))*(Ini[k1,k2]-2*Ini[k1,k2]+Ini[k1-1,k2]+Ini[k1,k2-1]-2*Ini[k1,k2]+Ini[k1,k2+1]);
										elseif I>2
											#println("SBAM");
											#println(size((Ini[end-1])));
											S[k1,k2]=(Delta^2/(h^2))*(Ini[end][k1+1,k2]-2*Ini[end][k1,k2]+Ini[end][k1-1,k2]+Ini[end][k1,k2-1]-2*Ini[end][k1,k2]+Ini[end][k1,k2+1])+2*Ini[end][k1,k2]-Ini[end-1][k1,k2];
										end
									end
								end
								if I > 1 
									XSLength = r;
									YSLength = r;
									for m in 1:r
										k2 = r*j-r+m
										S[(r*i)-r,k2]=Gamma(X2[(r*j)-r],Y2[k2])
									end
									for m in 1:r
										k1 = r*i-r+m
										S[k1,(r*j)-r+YSLength]=Gamma(X2[k1],Y2[(r*j)-r+YSLength]) 	
									end
									
								end
							elseif Z[i+1,j+1] == 1
							#   o....*
							#   |    :
							#   o....*
								println("Solving with 1 border, green O");
								H = X[2]-X[1];
								h = X2[2]-X2[1];
                                                                plot(X2[r*i-r]+(H/2),Y2[r*j-r]-(H/2),"go");
								for m in 1:r
									for n in 1:r
										k1 = (r*i)-r+m;
										k2 = (r*j)-r+n;
										if I == 1
											Ini[k1,k2] = f(X2[k1],Y2[k2]);	
											Ini2[k1,k2] = df(X2[k1],Y2[k2]);
										elseif I == 2
											S[k1,k2]=f(X2[k1],Y2[k2])+Delta*df(X2[k1],Y2[k2])+(Delta^2/(h^2))*(Ini[k1+1,k2]-2*Ini[k1,k2]+Ini[k1-1,k2]+Ini[k1,k2-1]-2*Ini[k1,k2]+Ini[k1,k2+1]);
										elseif I>2
											#println("SBAM");
											#println(size((Ini[end-1])));
											S[k1,k2]=(Delta^2/(h^2))*(Ini[end][k1+1,k2]-2*Ini[end][k1,k2]+Ini[end][k1-1,k2]+Ini[end][k1,k2-1]-2*Ini[end][k1,k2]+Ini[end][k1,k2+1])+2*Ini[end][k1,k2]-Ini[end-1][k1,k2];
										end
									end
								end
								if I > 1 
									XSLength = r;
									YSLength = r;	
									for m in 1:r
										k2=r*j-r+m;
										S[(r*j)-r,k2]=Gamma(X2[(r*j)-r],Y2[k2])
									end


								end
							end

						end
					elseif Z[i+1,j]==2
						
						if Z[i,j+1] == 1
							if Z[i+1,j+1]==1
							#   o____o
							#   :    :
							#   *....*
								println("Solving with 1 border, green O");
								H = X[2]-X[1];
                                                                h = X2[2]-X2[1];
                                                                plot(X2[r*i-r]+(H/2),Y2[r*j-r]-(H/2),"go");

								for m in 1:r
									for n in 1:r
										k1 = (r*i)-r+m;
										k2 = (r*j)-r+n;
										if I == 1
											Ini[k1,k2] = f(X2[k1],Y2[k2]);	
											Ini2[k1,k2] = df(X2[k1],Y2[k2]);
										elseif I == 2
											S[k1,k2]=f(X2[k1],Y2[k2])+Delta*df(X2[k1],Y2[k2])+(Delta^2/(h^2))*(Ini[k1+1,k2]-2*Ini[k1,k2]+Ini[k1-1,k2]+Ini[k1,k2-1]-2*Ini[k1,k2]+Ini[k1,k2+1]);
										elseif I>2
											#println("SBAM");
											#println(size((Ini[end-1])));
											S[k1,k2]=(Delta^2/(h^2))*(Ini[end][k1+1,k2]-2*Ini[end][k1,k2]+Ini[end][k2-1,k2+n]+Ini[end][k1,k2-1]-2*Ini[end][k1,k2]+Ini[end][k1,k2+1])+2*Ini[end][k1,k2]-Ini[end-1][k1,k2];
										end
									end
								end
								if I > 1 
									XSLength = r;
									YSLength = r;
									for m in 1:r
										k1=r*i-r+m;
										S[k1,(r*j)-r]=Gamma(X2[k1],Y2[(r*j)-r+1]) 	
										plot(X2[(r*i)-r+m],Y2[(r*j-r)],"y*");
									end
									
								end
							elseif Z[i+1,j+1] == 2
							#   o____o
							#   :    |
							#   *....o
								println("Solving with 2 border, blue O");
								H = X[2]-X[1];
                                                                h = X2[2]-X2[1];
                                                                plot(X2[r*i-r]+(H/2),Y2[r*j-r]-(H/2),"go");
								for m in 1:r
									for n in 1:r
										k1 = (r*i)-r+m;
										k2 = (r*j)-r+n;
										if I == 1
											Ini[k1,k2] = f(X2[k1],Y2[k2]);	
											Ini2[k1,k2] = df(X2[k1],Y2[k2]);
										elseif I == 2
											S[k1,k2]=f(X2[k1],Y2[k2])+Delta*df(X2[k1],Y2[k2])+(Delta^2/(h^2))*(Ini[k1+1,k2]-2*Ini[k1,k2]+Ini[k1-1,k2]+Ini[k1,k2-1]-2*Ini[k1,k2]+Ini[k1,k2+1]);

										elseif I>2
											#println("SBAM");
											#println(size((Ini[end-1])));
											S[k1,k2]=(Delta^2/(h^2))*(Ini[end][k1+1,k2]-2*Ini[end][k1,k2]+Ini[end][k1-1,k2]+Ini[end][k1,k2-1]-2*Ini[end][k1,k2]+Ini[end][k1,k2+1])+2*Ini[end][k1,k2]-Ini[end-1][k1,k2];
										end
									end
								end
								if I > 1 
									XSLength = r;
									YSLength = r;
									for m in 1:r
										k1=r*i-r+m;
										S[k1,(r*j)-r]=Gamma(X2[k1],Y2[(r*j)-1]) 	
									end
									for m in 1:r
										k2=r*j-r+m;
										S[(r*i)-r+XSLength,k2]=Gamma(X2[(r*i)-r+XSLength],Y2[k2]) 	
									end
									
								end

							end
						elseif Z[i,j+1] == 2
							if Z[i+1,j+1] == 2
							#   o____o
							#   |    |
							#   o____o
								println("Solving on full square good luck finding it...");
							end
							if Z[i+1,j+1] == 1
							#   o----o
							#   |    :
							#   o....*
								println("Solving with 2 border, blue O");
								H = X[2]-X[1];
								h = X2[2]-X2[1];                            	
                                                                plot(X2[r*i-r]+(H/2),Y2[r*j-r]-(H/2),"go");
								for m in 1:r
									for n in 1:r
										k1 = (r*i)-r+m;
										k2 = (r*j)-r+n;
										if I == 1
											Ini[k1,k2] = f(X2[k1],Y2[k2]);	
											Ini2[k1,k2] = df(X2[k1],Y2[k2]);
										elseif I == 2
											S[k1,k2]=f(X2[k1],Y2[k2])+Delta*df(X2[k1],Y2[k2])+(Delta^2/(h^2))*(Ini[k1+1,k2]-2*Ini[k1,k2]+Ini[k1-1,k2]+Ini[k1,k2-1]-2*Ini[k1,k2]+Ini[k1,k2+1]);

										elseif I>2
											#println("SBAM");
											#println(size((Ini[end-1])));
											S[k1,k2]=(Delta^2/h^2)*(Ini[end][k1+1,k2]-2*Ini[end][k1,k2]+Ini[end][k1-1,k2]+Ini[end][k1,k2-1]-2*Ini[end][k1,k2]+Ini[end][k1,k2+1])+2*Ini[end][k1,k2]-Ini[end-1][k1,k2];
										end
									end
								end
								if I > 1 
									XSLength = r;
									YSLength = r;
									for m in 1:r
										k1 = (r*i)-r+m;
										S[k1,(r*j)-r]=Gamma(X2[k1],Y2[(r*j)-r]) 	
									end
									for m in 1:r
										k2=(r*j)-r+m;
										S[(r*i)-r,k2]=Gamma(X2[(r*i)-r],Y2[k2]) 	
									end
									
								end

							end
						end

					end
				end

			end
		end
	
	if I == 2
		tmp = [];
		tmp = push!(tmp,Ini);
		tmp = push!(tmp,S);
		Ini = tmp;	
	elseif I > 2
		Ini = push!(Ini,S);
	end
	close()
	if I == 1
		matshow(Ini,origin="upper",extent=[-3,3,-3,3],vmax=2.5,vmin=-0.5);
	elseif I>1
		matshow(Ini[end],origin="upper",extent=[-3,3,-3,3],vmax=2.5,vmin=-0.5);
	end
	colorbar();
	title("Initial Data on the Mesh, with boder");
	for i in 1:size(PN)[1]
		plot(PC[i][1],PC[i][2],"g*");
	end	
	savefig(string("Snap/Snapshot",I,".png"));
	
	
	println("[INFO]: size S ",size(S)," r: ",r," size Z:", size(Z)," proportion: ",size(S)[1]/size(Z)[1]);
	println("[INFO]: X: ", X[20]," X2: ",X2[r*20-r])
	println("[INFO]: X: ", X[30]," X2: ",X2[r*30-r])
	#figure();
	#println(size(spacer(X,10))[1]);
	#plot(X,ones(N,1),"ro");
	#plot(spacer(X,10),ones(N*r-r,1),"g*");
	#msg = readline()
	#if msg == "Q" return; end
	sleep(0.1);
	close();
end
end
