function FDirichlet(q,f,D,I,N)
#This function solve Dirichlet differential problem of the type:
# -u''(x)+q(x)u(x)=f(x);
#where q(x) is a postive function
#N is the dimension of the gird.
	h=(D[2]-D[1])/(N+1); #Computing step size
	x=D[1]:h:D[2];#computing the x vector, where we evaluate f and h;
	d1=(-1)*ones(1,N-1);
	d2=[];
	for i in 2:length(x)-1
		d2=push!(d2,2+q(x[i])*(h^2));
	end
	d2=convert(Array{Float64}, d2);
	A=spdiagm(-1=>d1[1,:],0=>d2,1=>d1[1,:])
	A=(1/(h^2))*A;
	b=[];
	for i in 2:length(x)-1
		b=push!(b,f(x[i]));
	end
	b[1]=b[1]+(I[1]/(h^2));
	b[length(b)]=b[length(b)]+(I[2]/(h^2));
	u=A\b;
	y = [I[1]];
	y=append!(y,u);
	y=push!(y,I[2]);
	return y, x, h, A;
end
function FDNeumann(q,f,D,I,N)
#This function solve Neumann differential problem of the type:
# -u''(x)+q(x)u(x)=f(x);
#where q(x) is a postive function
#N is the dimension of the gird.
	h=(D[2]-D[1])/(N+1); #Computing step size
	x=D[1]:h:D[2];#computing the x vector, where we evaluate f and h;
	d1=(-1)*(1/(h^2))*ones(1,N+1);
	d2=[(1/(h^2))+q(x[1])/2];
	for i in 2:length(x)-1
		d2=push!(d2,(2/(h^2))+q(x[i]));
	end
	d2=push!(d2,(1/(h^2))+q(x[N+1])/2);
	d2=convert(Array{Float64}, d2);
	A=spdiagm(-1=>d1[1,:],0=>d2,1=>d1[1,:]);
	b=[((f(x[1])/2)-(I[1]/h))];
	for i in 2:length(x)-1
		b=push!(b,f(x[i]));
	end
	b=push!(b,(I[2]/h)+f(x[N+1])/2)
	y=A\b;
	return y, x, h, A;
end
function FDNeumannSlow(q,f,D,I,N)
#This function solve Neumann differential problem of the type:
# -u''(x)+q(x)u(x)=f(x);
#where q(x) is a postive function
#N is the dimension of the gird.
	h=(D[2]-D[1])/(N+1); #Computing step size
	x=D[1]:h:D[2];#computing the x vector, where we evaluate f and h;
	d1=(-1)*(1/(h^2))*ones(1,N+1);
	d2=[(1/(h^2))]
	for i in 2:length(x)-1
		d2=push!(d2,(2/(h^2))+q(x[i]));
	end
	d2=push!(d2,(1/(h^2)));
	d2=convert(Array{Float64}, d2);
	A=spdiagm(-1=>d1[1,:],0=>d2,1=>d1[1,:]);
	b=[(-(I[1]/h))];
	for i in 2:length(x)-1
		b=push!(b,f(x[i]));
	end
	b=push!(b,(I[2]/h))
	y=A\b;
	return y, x, h, A;
end
function tridiagsolverAmbro(a,b,c,y)
	n = length(a);
	r = c;
	u = zeros(1,n);
	l = zeros(1,n-1);
	u[1]=a[1];
	for i in 2:n
        	l[i-1] = b[i-1] / u[i-1];
        	u[i] = a[i] - l[i-1] * r[i-1];
    	end
	z = zeros(1,n);
	z[1] = y[1];
	for j in 2:n
        	z[j] = y[j] - l[j-1] * z[j-1];
    	end
	x = zeros(n, 1);
	x[n] = z[n] / u[n];
	for k in 1:n-1
       		x[n-k] = (z[n-k] - r[n-k] * x[n-k+1]) / u[n-k];
    	end
	return x;
end
function FDNonLinearDirchletAutonomus(f,df,I,D,N,u,tol=0.0001,it=100)
	h=(I[2]-I[1])/(N+1);
	x=I[1]:h:I[2];
	for k in 1:it
		G=[]
		DG=[];
		#Computing vector G
		G=[(D[1]-2*u[1]+u[2])/(h^2)+f(u[1])];
		for j in 2:N-1
			G=push!(G,(u[j-1]-2*u[j]+u[j+1])/(h^2)+f(u[j]));
		end
		G=push!(G,(u[N-1]-2*u[N]+D[2])/(h^2)+f(u[N]));
		d1=(1/(h^2))*ones(1,N-1);
		d2=[];
		for i in 1:N
			d2=push!(d2,(-2/(h^2))+df(u[i]));
		end
		w = tridiagsolverAmbro(d2,d1,d1,G);
		u=u-w';
		if (norm(G,2)<tol)
			y = [D[1]];
			y = append!(y,u);
			y = push!(y,D[2]);
	return x, y;
		end

	end
	y = [D[1]];
	y = append!(y,u);
	y = push!(y,D[2]);
	return x, y;

end
function FEM1H(q,f,D,N) #FiniteElementMethd1inearHomogeneus
	h=(D[2]-D[1])/(N+1);
	x=D[1]:h:D[2];
	println(length(x));
	#Costruisco il vettore Q, che è la funzione q valutata nei punti medi.
	Q=[];
	for i in 2:N+2
		Q=push!(Q,q((x[i-1]+x[i])/2));
	end
	println(length(Q));
	#Costruiamo le varie diagonali
	d2=[];
	for i in 1:N
		d2=push!(d2,(2/h)+(Q[i]*h+Q[i+1]*h)/3);
	end
	d1=[];
	for i in 2:N
		d1=push!(d1,-(1/h)+(Q[i]*h)/6);
	end
	println(length(d2));
	println(length(d1));
	ff=[];
	for i in 2:N+2
		ff=push!(ff,f((x[i-1]+x[i])/2));
	end
	F = [];
	for i in 1:N
		F=push!(F,(h*ff[i]+h*ff[i+1])/2);
	end
	y=tridiagsolverAmbro(d2,d1,d1,F);
	s = [0.0];
	s=append!(s,y);
	s=push!(s,0.0);
	return s, x;
end
function FEM1(q,f,D,I,N)
	#Costruisco funzione di lifting l
	l(x) =((I[2]-I[1])/(D[2]-D[1]))*x-((I[2]-I[1])/(D[2]-D[1]))*D[1]+I[1];
	#Risolvo il problema: -u''+qu=f+u''-qu
	F(x)=f(x)-q(x)*l(x);
	y, x =FEM1H(q,F,D,N);
	println("length FEM1H:",length(y));
	s=[];
	println(x);
	println(y);
	for i in 1:length(x)
		s=push!(s,y[i]+l(x[i]));
	end
	print(s);
	return s, x;
end
function FDTimeStepHeat(f,D,Ini,N,r,t;opt="Euler")
	#We will try here to solve the differential equation of heat diffusion.
	if opt=="Euler"
		h=(D[2]-D[1])/(N+1); #Computing step size
		x=D[1]:h:D[2];#computing the x vector, where we evaluate f and h;
		U=[];
		for i in 2:length(x)-1
			U=push!(U,f(x[i]));
		end
		U = convert(Array{Float64}, U);
		d1=(-1)*ones(1,N-1);
		d2=[];
		for i in 2:length(x)-1
			d2=push!(d2,2);
		end
		d2=convert(Array{Float64}, d2);
		A=spdiagm(-1=>d1[1,:],0=>d2,1=>d1[1,:]);
		A=(1/(h^2))*A;
		#println(size(U));
		#println(size(x));
		T = 0;
		Δt = h^r;
		while true
			T = T+Δt;
			U=(sparse(I,N,N)-(h^(r+2))*A)*U;
			Y = [f(D[1])];
			Y = append!(Y,U);
			Y = push!(Y,f(D[2]));
			if (t<T)
				return Y, h
			end
		end
	end
end
function LeapfrogWave(f,g,N,D,Ini,r,t,v=1)
	#We will try here to solve the differential equation of wave diffusion.
	#N=100;
	#D=[0,1];
	h=(D[2]-D[1])/(N+1); #Computing step size
	ΔT = r*h;
	x=D[1]:h:D[2];#computing the x vector, where we evaluate f and h;
	#f(x)=x*(1-x);
	#g(x)=0*x;

	V=[];
	U=[];
	for i in 2:length(x)-1
		U=push!(U,f(x[i]));
	end
	U = convert(Array{Float64}, U);
	V = push!(V,U);

	Y = [Ini[1]];
	Y = append!(Y,V[1]);
	Y = push!(Y,Ini[2]);
	U=[];
	for i in 2:length(x)-1
		U=push!(U,Y[i]+ΔT*g(x[i])+v^2*(ΔT^2/(2*h^2))*(Y[i-1]-2Y[i]+Y[i+1]));
	end
	U = convert(Array{Float64}, U);
	V = push!(V,U);

	d1=(-1)*ones(1,N-1);
	d2=[];
	for i in 2:length(x)-1
		d2=push!(d2,2);
	end
	d2=convert(Array{Float64}, d2);
	A=spdiagm(-1=>d1[1,:],0=>d2,1=>d1[1,:]);
	A=(v^2/(h^2))*A;
	Y = [Ini[1]];
	Y = append!(Y,V[1]);
	Y = push!(Y,Ini[2]);
	T=ΔT;
	Out=[];
	while true
		if length(t)==1
			T = T+ΔT;
			V=push!(V(2*sparse(I,N,N)-(ΔT)^2*A)*V[end]-V[end-1]);
			Y = [Ini[1]];
			Y = append!(Y,V[end]);
			Y = push!(Y,Ini[2]);
			if (t<T)
				return Y, h
			end
		elseif length(t)==2
			T = T+ΔT;
			V=push!(V(2*sparse(I,N,N)-(ΔT)^2*A)*V[end]-V[end-1]);
			Y = [Ini[1]];
			Y = append!(Y,V[end]);
			Y = push!(Y,Ini[2]);
			Out = push!(Out,Y)
			if (t[2]<T)
				return Out, h
			end
		end
	end
end
function LeapfrogWaveMedium(f,g,N,D,Ini,r,t,v)
	#We will try here to solve the differential equation of wave diffusion.
	#N=100;
	#D=[0,1];
	h=(D[2]-D[1])/(N+1); #Computing step size
	ΔT = r*h;
	x=D[1]:h:D[2];#computing the x vector, where we evaluate f and h;
	#f(x)=x*(1-x);
	#g(x)=0*x;

	V=[];
	U=[];
	for i in 2:length(x)-1
		U=push!(U,f(x[i]));
	end
	U = convert(Array{Float64}, U);
	V = push!(V,U);

	Y = [Ini[1]];
	Y = append!(Y,V[1]);
	Y = push!(Y,Ini[2]);
	U=[];
	for i in 2:length(x)-1
		U=push!(U,Y[i]+ΔT*g(x[i])+(v(x[i])^2)*(ΔT^2/(2*h^2))*(Y[i-1]-2Y[i]+Y[i+1]));
	end
	U = convert(Array{Float64}, U);
	V = push!(V,U);

	d1=[];
	for i in 2:length(x)-2
		d1=push!(d1,(-1)*v(x[i])^2);	
	end
	d1=convert(Array{Float64}, d1);
	
	d2=[];
	for i in 2:length(x)-1
		d2=push!(d2,2*v(x[i])^2);
	end
	d2=convert(Array{Float64}, d2);
	A=spdiagm(-1=>d1,0=>d2,1=>d1);
	#print(A);
	A=(1/(h^2))*A;
	Y = [Ini[1]];
	Y = append!(Y,V[1]);
	Y = push!(Y,Ini[2]);
	T=ΔT;
	Out=[];
	while true
		if length(t)==1
			T = T+ΔT;
			V=push!(V,(2*sparse(I,N,N)-(ΔT)^2*A)*V[end]-V[end-1]);
			Y = [Ini[1]];
			Y = append!(Y,V[end]);
			Y = push!(Y,Ini[2]);
			if (t<T)
				return Y, h
			end
		elseif length(t)==2
			T = T+ΔT;
			V=push!(V,(2*sparse(I,N,N)-(ΔT)^2*A)*V[end]-V[end-1]);
			Y = [Ini[1]];
			Y = append!(Y,V[end]);
			Y = push!(Y,Ini[2]);
			Out = push!(Out,Y)
			if (t[2]<T)
				return Out, h
			end
		end
	end
end
function NeumannLeapfrogWave(f,g,N,D,Ini,r,t)
	h=(D[2]-D[1])/(N+1); #Computing step size
	ΔT = r*h;
	x=D[1]:h:D[2];#computing the x vector, where we evaluate f and h;

	V=[];
	U=[];
	for i in 2:length(x)-1
		U=push!(U,f(x[i]));
	end
	U = convert(Array{Float64}, U);
	V = push!(V,U);

	Y = [Ini[1]];
	Y = append!(Y,V[1]);
	Y = push!(Y,Ini[2]);
	U=[];
	for i in 2:length(x)-1
		U=push!(U,Y[i]+ΔT*g(x[i])+(ΔT^2/(2*h^2))*(Y[i-1]-2Y[i]+Y[i+1]));
	end
	U = convert(Array{Float64}, U);
	V = push!(V,U);

	d1=(-1)*ones(1,N-1);
	d2=[];
	for i in 2:length(x)-1
		d2=push!(d2,2);
	end
	d2=convert(Array{Float64}, d2);
	A=spdiagm(-1=>d1[1,:],0=>d2,1=>d1[1,:]);
	A=(1/(h^2))*A;
	Y = [Ini[1]];
	Y = append!(Y,V[1]);
	Y = push!(Y,Ini[2]);
	T=0;
	while true
		T = T+ΔT;
		V=push!(V,(2*sparse(I,N,N)-(ΔT)^2*A)*V[end]-V[end-1]);
		Y = [Ini[1]];
		Y = append!(Y,V[end]);
		Y = push!(Y,Ini[2]);
		if (t<T)
			return Y, h
		end
	end
end
function ForcedLeapfrogWave(f,g,F,N,D,Ini,r,t)
	#We will try here to solve the differential equation of wave diffusion.
	#N=100;
	#D=[0,1];
	h=(D[2]-D[1])/(N+1); #Computing step size
	ΔT = r*h;
	x=D[1]:h:D[2];#computing the x vector, where we evaluate f and h;
	#f(x)=x*(1-x);
	#g(x)=0*x;

	V=[];
	U=[];
	for i in 2:length(x)-1
		U=push!(U,f(x[i]));
	end
	U = convert(Array{Float64}, U);
	V = push!(V,U);

	Y = [Ini[1]];
	Y = append!(Y,V[1]);
	Y = push!(Y,Ini[2]);
	U=[];
	for i in 2:length(x)-1
		U=push!(U,Y[i]+ΔT*g(x[i])+(ΔT^2/(2*h^2))*(Y[i-1]-2Y[i]+Y[i+1]));
	end
	U = convert(Array{Float64}, U);
	V = push!(V,U);

	d1=(-1)*ones(1,N-1);
	d2=[];
	for i in 2:length(x)-1
		d2=push!(d2,2);
	end
	d2=convert(Array{Float64}, d2);
	A=spdiagm(-1=>d1[1,:],0=>d2,1=>d1[1,:]);
	A=(1/(h^2))*A;
	Y = [Ini[1]];
	Y = append!(Y,V[1]);
	Y = push!(Y,Ini[2]);
	T=0;
	while true
		T = T+ΔT;
		V=push!(V,(2*sparse(I,N,N)-(ΔT)^2*A)*V[end]-V[end-1]+F.(x,T)[2:end-1]);
		Y = [Ini[1]];
		Y = append!(Y,V[end]);
		Y = push!(Y,Ini[2]);
		if (t<T)
			return Y, h
		end
	end
end
function HLeapFrogWave2D(f,∂f,N,D,r,t)
	Δx = (D[2]-D[1])/(N+1);
	Δt = Δx*r;
	Y=[];
	φ = ones(size(D)[1],size(D)[1])
	for i in 1:size(D)[1]
		for j in 1:size(D)[1]
			φ[i,j] = f(D[i],D[j]);
		end
	end
	φ2 = ones(size(D)[1],size(D)[1])
	for i in 1:size(D)[1]
		for j in 1:size(D)[1]
			if i == 1 || j == 1
				φ2[i,j]=0.0;
			elseif i==size(D)[1] || j==size(D)[1]
				φ2[i,j]=0.0;
			else
				φ2[i,j] = f(D[i],D[j])+ Δt * ∂f(D[i],D[j])+((Δt^2)/(Δx^2))*(φ[i+1,j]-2*φ[i,j]+φ[i-1,j]+φ[i,j+1]-2*φ[i,j]+φ[i,j-1]);
			end
		end
	end
	Y = push!(Y,φ);
	Y = push!(Y,φ2);
	T = 0;
	while true
		T = T+Δt;
		Y2 = ones(size(D)[1],size(D)[1]);
		for i in 1:size(D)[1]
			for j in 1:size(D)[1]
				if i == 1 || j == 1
			  		Y2[i,j]=0.0;
			  	elseif i==size(D)[1] || j==size(D)[1]
			  		Y2[i,j]=0.0;
			  	else
			  		Y2[i,j] = (Δt^2/Δx^2)*(Y[end][i+1,j]-2*Y[end][i,j]+Y[end][i-1,j]+Y[end][i,j+1]-2*Y[end][i,j]+Y[end][i,j-1])+2*Y[end][i,j]-Y[end-1][i,j];
			  	end
			end
		end
		Y = push!(Y,Y2);
		if (t<T)
			return Y2, Δx;
		end
	end
end

function MillerWave(f,g,N,D,Ini,r,t,v=1,R=1)

	h=(D[2]-D[1])/(N+1); #Computing step size
	Δt = r*h;
	x=D[1]:h:D[2];
	U=[];
	R=r;
	#Aggiungiamo dato iniziale
	U= push!(U,f.(x));
	V = [Ini[1]]
	for i in 2:length(U[1])-1
		V=push!(V,U[1][i]+Δt*g(x[i])+0.5*v^2*R^2*(U[1][i-1]-2*U[1][i]+U[1][i+1]))
	end
	V = push!(V,Ini[2]);
	U=push!(U,V);

	U[1] = deleteat!(U[1],1);
	U[1] = deleteat!(U[1],length(U[1]));
	U[2] = deleteat!(U[2],1);
	U[2] = deleteat!(U[2],length(U[2]));

	θ=(1+(2/(v^2*r^2)));
	γ=(1-(2/(v^2*r^2)))

	d1=ones(1,N-1);
	d1=convert(Array{Float64}, d1);
	d2=2*ones(1,N);
	d2=convert(Array{Float64}, d2);

	A=spdiagm(-1=>(-1)*d1[1,:],0=>θ*d2[1,:],1=>(-1)*d1[1,:]);
	B=spdiagm(-1=>d1[1,:],0=>-γ*d2[1,:],1=>d1[1,:]);
	C=spdiagm(-1=>d1[1,:],0=>-θ*d2[1,:],1=>d1[1,:]);

	T = Δt;
	Out = [];
	while true
		if length(t)==1
			T=T+Δt;
			b = 2*B*U[end]+C*U[end-1];
			U = push!(U,A\b);
			Y = [Ini[1]];
			Y = append!(Y,U[end]);
			Y = push!(Y,Ini[2]);
			if (t< T)
				return Y,h;
			end
		elseif length(t)==2
			T=T+Δt;
			b = 2*B*U[end]+C*U[end-1];
			U = push!(U,A\b);
			Y = [Ini[1]];
			Y = append!(Y,U[end]);
			Y = push!(Y,Ini[2]);
			Out=push!(Out,Y)
			if (t[2]< T)
				return Out,h;
			end
		end
	end
end
function MillerZWave(f,g,N,D,Ini,r,t,v=1,R=1)

	h=(D[2]-D[1])/(N+1); #Computing step size
	Δt = r*h;
	x=D[1]:h:D[2];
	U=[];
	R=r;
	#Aggiungiamo dato iniziale
	U= push!(U,f.(x));
	V = [Ini[1]]
	for i in 2:length(U[1])-1
		V=push!(V,U[1][i]+Δt*g(x[i])+0.5*R^2*(U[1][i-1]-2*U[1][i]+U[1][i+1]))
	end
	V = push!(V,Ini[2]);
	U=push!(U,V);

	U[1] = deleteat!(U[1],1);
	U[1] = deleteat!(U[1],length(U[1]));
	U[2] = deleteat!(U[2],1);
	U[2] = deleteat!(U[2],length(U[2]));

	θ=(1+(2/(v^2*r^2)));
	γ=(1-(2/(v^2*r^2)))

	d1=ones(1,N-1);
	d1=convert(Array{Float64}, d1);
	d2=2*ones(1,N);
	d2=convert(Array{Float64}, d2);

	A=spdiagm(-1=>(-1)*d1[1,:],0=>θ*d2[1,:],1=>(-1)*d1[1,:]);
	B=spdiagm(-1=>d1[1,:],0=>-γ*d2[1,:],1=>d1[1,:]);
	C = 2*inv(Matrix(A))*B;

	T = Δt;
	Out = [];
	while true
		if length(t)==1
			T=T+Δt;
			Δ = C*U[end];
			U = push!(U,Δ-U[end-1]);
			Y = [Ini[1]];
			Y = append!(Y,U[end]);
			Y = push!(Y,Ini[2]);
			if (t< T)
				return Y,h;
			end
		elseif length(t)==2
			T=T+Δt;
			Δ = C*U[end];
			U = push!(U,Δ-U[end-1]);
			Y = [Ini[1]];
			Y = append!(Y,U[end]);
			Y = push!(Y,Ini[2]);
			Out=push!(Out,Y)
			if (t[2]< T)
				return Out,h;
			end
		end
	end
end
function NewmarkWave(f,g,N,D,Ini,r,t,v=1)

	h=(D[2]-D[1])/(N+1); #Computing step size
	Δt = r*h;
	x=D[1]:h:D[2];
	U=[];
	V=[];
	A=[];

	d1=ones(1,N-1);
	d1=convert(Array{Float64}, d1);
	d2=-2*ones(1,N);
	d2=convert(Array{Float64}, d2);

	DFC=spdiagm(-1=>(1)*d1[1,:],0=>d2[1,:],1=>(1)*d1[1,:]);
	DFC=((-v^2)/(h^2))*DFC;
	#Aggiungiamo dato iniziale
	U = push!(U,f.(x));
	V = push!(V,g.(x));

	U[1] = deleteat!(U[1],1);
	U[1] = deleteat!(U[1],length(U[1]));
	V[1] = deleteat!(V[1],1);
	V[1] = deleteat!(V[1],length(V[1]));

	A = push!(A,-DFC*U[end]);

	β=0.25;
	γ=0.5;

	T=0;
	Out=[];
	while true
		if length(t)==1
			T=T+Δt;
			A = push!(A,(sparse(I,N,N)+β*(Δt^2)*DFC)\((-DFC)*(U[end]+Δt*V[end]+Δt^2*(0.5-β)*A[end])))
			U = push!(U,U[end]+Δt*V[end]+Δt^2*((0.5-β)*A[end-1]+β*A[end]));
			V = push!(V,V[end]+Δt*((1-γ)*A[end-1] + γ*A[end]));
			Y = [Ini[1]];
			Y = append!(Y,U[end]);
			Y = push!(Y,Ini[2]);
			if (t< T)
				return Y,h;
			end
		elseif length(t)==2
			T=T+Δt;
			A = push!(A,(sparse(I,N,N)+β*(Δt^2)*DFC)\((-DFC)*(U[end]+Δt*V[end]+Δt^2*(0.5-β)*A[end])))
			U = push!(U,U[end]+Δt*V[end]+Δt^2*((0.5-β)*A[end-1]+β*A[end]));
			V = push!(V,V[end]+Δt*((1-γ)*A[end-1] + γ*A[end]));
			Y = [Ini[1]];
			Y = append!(Y,U[end]);
			Y = push!(Y,Ini[2]);
			Out=push!(Out,Y)
			if (t[2]< T)
				return Out,h;
			end
		end
	end
end
