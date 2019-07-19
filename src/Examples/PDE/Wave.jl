function waveExample(v=1,Auto=true)
	#We will try here to solve the differential equation of wave diffusion.
	N=100;
	D=[0.0,1.0];
	Ini=[0.0,0.0];
	f(x)=sin(π*x);
	g(x)=0*x;
	e(x,t)=sin(π*x)*cos(π*t*v)
	h=(D[2]-D[1])/(N+1);
	r=0.5;
	t = r*h;
	err = [];
	while 0 < t <π/2
		if Auto
			sleep(0.1);
		else
			readline();
		end
		PyPlot.clf();
		ax = gca()
		ax[:set_ylim]([-1.3,1.3])
		Y,h = LeapfrogWave(f,g,N,D,Ini,r,t,v);
		Δt = r*h;
		t=t+Δt;
		x=D[1]:h:D[2];
		E = e.(x,t);
		err = push!(err,norm(Y-E,Inf));
		plot(x,Y)
		plot(x,E)
		plot([0.0,1.0],[1.0,1.0],"r-");
		title(string("Time: ",t," r:",r));
	end
	return(err);
end
function TsunamiWaveExample(v=1,Auto=true)
	#We will try here to solve the differential equation of wave diffusion.
	N=100;
	D=[-2.0,2.0];
	Ini=[0.0,0.0];
	μ=0.4; σ=0.1;

	f(x)=exp((-(x-μ)^2)/(2*σ^2));
	g(x)=0;
	#e(x,t)=sin(π*x)*cos(π*t*v)
	t = 0;
	Δt = 0.1;
	while t < 4
		if Auto
			sleep(0.5);
		else
			readline();
		end
		PyPlot.clf();
		ax = gca()
		ax[:set_ylim]([-1.3,1.3])
		Y,h = LeapfrogWave(f,g,N,D,Ini,0.35,t,v);
		x=D[1]:h:D[2];
		#E = e.(x,t);
		t=t+Δt;
		plot(x,Y)
		title(string("Time: ",t));
	end

end
function WaveMediumExample(Auto=true)
	#We will try here to solve the differential equation of wave diffusion.
	v(x) = if (0 ≤ x ≤  0.5) return 0.5  else return 0.8 end
	N=100;
	D=[0.0,1.0];
	Ini=[0.0,0.0];
	f(x)=exp((-(x-0.2)^2)/(2*0.001));
	g(x)=0*x;
	h=(D[2]-D[1])/(N+1);
	r=0.2;
	t = r*h;
	err = [];
	while 0 < t < 2
		if Auto
			sleep(0.1);
		else
			readline();
		end
		PyPlot.clf();
		ax = gca()
		ax[:set_ylim]([-1.3,1.3])
		Y,h = LeapfrogWaveMedium(f,g,N,D,Ini,r,t,v);
		Δt = r*h;
		t=t+Δt;
		x=D[1]:h:D[2];
		plot(x,Y)
		#plot(x,v.(x))
		plot([0.5,0.5],[-1,1],"r--");
		#plot([0.0,1.0],[1.0,1.0],"r-");
		title(string("Time: ",t," r:",r));
	end
	return(err);
end
function MillerWaveExample(v=1,Auto=true)
	N=100;
	D=[0.0,1.0];
	Ini=[0.0,0.0];
	f(x)=sin(π*x);
	e(x,t)=sin(π*x)*cos(π*t*v)
	g(x)=0*x;
	r = 1.0;
	R = 1;

	h=(D[2]-D[1])/(N+1); #Computing step size
	Δt = r*h;
	x=D[1]:h:D[2];
	U=[];
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
	C=spdiagm(-1=>d1[1,:],0=>-θ*d2[1,:],1=>d1[1,:]);

	for t in 2*Δt:Δt:4
		b = 2*B*U[end]+C*U[end-1];
		U = push!(U,A\b);
		sleep(0.005);
		PyPlot.clf();
		ax = gca()
		ax[:set_ylim]([-1.2,1.2])
		Y = [Ini[1]];
		Y = append!(Y,U[end]);
		Y = push!(Y,Ini[2]);
		x=D[1]:h:D[2];
		E = e.(x,t);
		t=t+Δt;
		plot(x,Y)
		plot(x,E)
		plot([0.0,1.0],[1.0,1.0],"r-");
		title(string("Time: ",t));
	end
end
function MillerZWaveExample(v=1,Auto=true)
	N=100;
	D=[0.0,1.0];
	Ini=[0.0,0.0];
	f(x)=sin(π*x);
	e(x,t)=sin(π*x)*cos(π*t*v)
	g(x)=0*x;
	r = 1.0;
	R = 1;

	h=(D[2]-D[1])/(N+1); #Computing step size
	Δt = r*h;
	x=D[1]:h:D[2];
	U=[];
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

	for t in 2*Δt:Δt:4
		Δ = C*U[end];
		U = push!(U,Δ-U[end-1]);
		sleep(0.005);
		PyPlot.clf();
		ax = gca()
		ax[:set_ylim]([-1.2,1.2])
		Y = [Ini[1]];
		Y = append!(Y,U[end]);
		Y = push!(Y,Ini[2]);
		x=D[1]:h:D[2];
		E = e.(x,t);
		t=t+Δt;
		plot(x,Y)
		plot(x,E)
		plot([0.0,1.0],[1.0,1.0],"r-");
		title(string("Time: ",t));
	end
end
function NewmarkWaveExample(Auto=true)
	N=100;
	D=[0.0,1.0];
	Ini=[0.0,0.0];
	f(x)=x*(1-x);
	g(x)=0*x;
	r = 1;


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
	DFC=((-1)/(h^2))*DFC;
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
	for t in 0:Δt:4
		A = push!(A,(sparse(I,N,N)+β*(Δt^2)*DFC)\((-DFC)*(U[end]+Δt*V[end]+Δt^2*(0.5-β)*A[end])))
		U = push!(U,U[end]+Δt*V[end]+Δt^2*((0.5-β)*A[end-1]+β*A[end]));
		V = push!(V,V[end]+Δt*((1-γ)*A[end-1] + γ*A[end]));
		sleep(0.05);
		PyPlot.clf();
		ax = gca()
		ax[:set_ylim]([-0.3,0.3])
		Y = [Ini[1]];
		Y = append!(Y,U[end]);
		Y = push!(Y,Ini[2]);
		plot(x,Y);
		title(t);
		plot([0.0,1.0],[0.25,0.25],"r-");
	end
end
function ForcedWaveExample(Auto=true)
	#We will try here to solve the differential equation of wave diffusion.
	N=100;
	D=[0.0,1];
	Ini=[0.0,0.0];
	σ=sqrt(2);
	F(x,t)=0.001*(x-0.5)*exp((-(t-10)^2)/(2*σ^2));
	f(x)=x*(x-1);
	g(x)=0*x;
	t = 0;
	Δt = 0.1;
	while true
		if Auto
			sleep(0.1);
		else
			readline();
		end
		PyPlot.clf();
		ax = gca()
		ax[:set_ylim]([-0.3,0.3])
		Y,h = ForcedLeapfrogWave(f,g,F,N,D,Ini,1,t);
		x=D[1]:h:D[2];
		t=t+Δt;
		plot(x,Y)
		title(string("Time: ",t));
	end
end
function waveExample2D(Auto=true)
	#f(x,y)=(x*(1-x))*(y*(1-y));
	μx=0; μy=1; σ=0.1;
	μx2=0; μy2=-1;
	f(x,y)=exp((-(x-μx)^2)/(2*σ^2)) * exp((-(y-μy)^2)/(2*σ^2))+exp((-(x-μx2)^2)/(2*σ^2)) * exp((-(y-μy2)^2)/(2*σ^2));
	∂f(x,y)=0;
	D=-4:0.05:4;
	N=100;
	t=0;
	Δt=0.001;
	while true
		if Auto
			sleep(0.01);
		else
			readline();
		end
		PyPlot.clf();
		ax = gca()
		Y,Δx = HLeapFrogWave2D(f,∂f,N,D,0.2,t);
		t=t+Δt;
		Zlimp = 0.1*ones(size(D)[1],size(D)[1]);
		Zlimm = -0.1*ones(size(D)[1],size(D)[1]);
		#plot_surface(D, D, Y, rstride=2,edgecolors="k", cstride=2,cmap=ColorMap("plasma"), alpha=0.8, linewidth=0.25);
		imshow(Y,vmin=-0.25,vmax=0.25);
		colorbar();
		#plot_surface(D, D, Zlimp, rstride=2,edgecolors="k", cstride=2,cmap=ColorMap("plasma"), alpha=0.8, linewidth=0.25);
		#plot_surface(D, D, Zlimm, rstride=2,edgecolors="k", cstride=2,cmap=ColorMap("plasma"), alpha=0.8, linewidth=0.25);
		title(string("Time: ",t));
		savefig(string("./Snap/t",t,".png"));
	end
end
function FDTDExample(Auto=true)
	#Here we simulate the right hand side of the wave graph (before node in strings).
	D=[0,2]; T=[0,10];
	N=100;
	h = (D[2]-D[1])/(N+1);
	r=1.00
	Δt = r*h;
	S = D[1]:(h/2):D[2];
	τ=1.5; σ=0.1;
	μ =1; ϵ=1;
	G(t)=exp((-(t-τ)^2)/(2*σ^2)); #Gaussian Pulse
	#G(t)=exp((-(t-τ)^2)/(2*σ^2))+exp((-(t-τ-4)^2)/(2*σ^2)); #Double Gaussian In Phase
	#G(t)=exp((-(t-τ)^2)/(2*σ^2))-exp((-(t-τ-4)^2)/(2*σ^2)); #Opposite Gaussian
	#G(t)=exp((-(t-τ)^2)/(2*σ^2))-exp((-(t-τ-1)^2)/(2*σ^2)); #Duble not in phase Gaussian

	X=[]; Y=[];
	for i in 1:length(S)
		if mod(i,2)==0
			Y = push!(Y,S[i]);
			#plot(i,S[i],"r*");
		else
			X = push!(X,S[i]);
			#plot(i,S[i],"b*");
		end
	end
	t=T[1]:Δt:T[2];
	figure();
	plot(t,G.(t));
	title("Gaussian Pulse");
	figure();
	E = zeros(length(X));
	H = zeros(length(Y));
	for dt in t
		if Auto
			sleep(0.005);
		else
			readline();
		end
		E[1]=G(dt);
		for i in 1:length(Y)
			H[i]=H[i]+(Δt/μ)*((E[i+1]-E[i])/h)
		end
		for i in 2:length(Y);
			E[i]=E[i]+(Δt/ϵ)*((H[i]-H[i-1])/h)
		end
		PyPlot.clf();
		ax=gca();
		ax[:set_ylim]([-2.3,2.3])
		title(string("Time: ",dt));
		plot(X,E);
		plot(Y,H);
	end
end
function CNFDTDExample(Auto=true)
    #Here we simulate the right hand side of the wave graph (before node in strings).
	D=[0,2]; T=[0,4];
	N=100;
	h = (D[2]-D[1])/(N-1);
	r=2;
	Δt = r*h;
	τ=0.4; σ=0.1;
    θ=0.5;
	G(t)=exp((-(t-τ)^2)/(2*σ^2)); #Gaussian Pulse

	x=D[1]:h:D[2];
	figure();

	E = zeros(length(x),1);
	H = zeros(length(x),1);

    #Costruiamo la matrice delle differenze finite
    d1=(-1)*ones(1,N-1);
    d2=[];
    d3=ones(1,N);
    for i in 2:length(x)-1
        d2=push!(d2,2.0);
    end
    d2=convert(Array{Float64}, d2);
    A=spdiagm(-1=>d1[1,:],0=>d2,1=>d1[1,:]);
    A=(1/(h^2))*A;
    B=spdiagm(-1=>d1[1,:],0=>d3[1,:]);
    B=(1/(h))*B;

    E = G.(x)
    #H = G.(x)
	for t in T[1]:Δt:T[2]
	 	if Auto
	 		sleep(0.05);
	 	else
	 		readline();
	 	end

	 	PyPlot.clf();
	 	ax=gca();
	 	ax[:set_ylim]([-2.3,2.3]);
	 	plot(x,E);
	 	plot(x,H);
        plot([0,2.0],[1,1],"r-");
        title(string("Time: ",t));

        #println(I)

        b=(H - Δt*B*E - Δt^2*θ*(1-θ)*A*H);
        H = (sparse(I,N,N)+Δt^2*θ^2*A)\b;

        b=(E - Δt*B*H - Δt^2*θ*(1-θ)*A*E);
        E = (sparse(I,N,N)+Δt^2*θ^2*A)\b;
	end
end
function FDTDMediumExample(Auto=true)
	D=[0,3]; T=[0,5];
	N=100;
	h = (D[2]-D[1])/(N+1);
	r=0.95;
	Δt = r*h;
	S = D[1]:(h/2):D[2];
	τ=1.5; σ=0.1;
	G(t)=exp((-(t-τ)^2)/(2*σ^2)); #Gaussian Pulse
	#G(t)=exp((-(t-τ)^2)/(2*σ^2))+exp((-(t-τ-4)^2)/(2*σ^2)); #Double Gaussian In Phase
	#G(t)=exp((-(t-τ)^2)/(2*σ^2))-exp((-(t-τ-4)^2)/(2*σ^2)); #Opposite Gaussian
	#G(t)=exp((-(t-τ)^2)/(2*σ^2))-exp((-(t-τ-1)^2)/(2*σ^2)); #Duble not in phase Gaussian
	function M(x)
		if x<1
			return 1
		else
			return 2
		end
	end
	function E(x)
		if x<1
			return 1
		else
			return 0.5
		end
	end
	X=[]; Y=[];
	for i in 1:length(S)
		if mod(i,2)==0
			Y = push!(Y,S[i]);
			#plot(i,S[i],"r*");
		else
			X = push!(X,S[i]);
			#plot(i,S[i],"b*");
		end
	end
	μ = M.(Y);
	ϵ = E.(Y);
	println(μ)
	t=T[1]:Δt:T[2];
	figure();
	plot(t,G.(t));
	title("Gaussian Pulse");
	figure();
	E = zeros(length(X));
	H = zeros(length(Y));
	for dt in t
		if Auto
			sleep(0.05);
		else
			readline();
		end
		E[1]=G(dt);
		for i in 1:length(Y)
			H[i]=H[i]-(Δt/μ[i])*((E[i+1]-E[i])/h)
		end
		for i in 2:length(Y)
			E[i]=E[i]-(Δt/ϵ[i])*((H[i]-H[i-1])/h)
		end
		println(E[end]);
		println(H[end])
		PyPlot.clf();
		title(string("Time: ",dt));
		ax=gca();
		ax[:set_ylim]([-2.3,2.3])
		plot([1,1],[-2.3,2.3],"r-");
		plot(X,E);
		plot(Y,H);
	end
end
