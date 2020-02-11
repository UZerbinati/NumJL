function LeapfrogBanchTimeWindow(r,T)
	D=[0.0,1.0];
	Ini=[0.0,0.0];
	f(x)=sin(π*x);
	g(x)=0*x;
	e(x,t)=sin(π*x)*cos(π*t)
	Energy = [];
	MaxIT = 1010;
	Error = zeros(Int64(MaxIT/10),5);
    i=1;
    for N in 10:10:MaxIT
        for j in 1:20
            τ = @elapsed begin
                Y,h = LeapfrogWave(f,g,N,D,Ini,r,T);
            end
            E = []
	    Amp = [];
            x=D[1]:h:D[2];
            t=r*h;
	    for k in 2:length(Y)-1
		y = Y[k];
		t=t+r*h;
		A = maximum(e.(x,t));
                E = push!(E,norm(y-e.(x,t),Inf))
		Amp = push!(Amp,abs(abs(A)-abs(maximum(y))));
		Energy = push!(Energy,abs(0.5*(π^2)-DirichletEnergy(Y[k-1],Y[k+1],h,r*h))); 
	    end
    	    x=D[1]:h:D[2];
            Error[i,1] = Error[i,1]+h;
            Error[i,2] = Error[i,2]+norm(E,Inf);
            Error[i,3] = Error[i,3]+τ
	    Error[i,4] = Error[i,4]+Amp[end];
	    Error[i,5] = Error[i,5]+Energy[end];
        end

        Error[i,1] = Error[i,1]/20;
        Error[i,2] = Error[i,2]/20;
        Error[i,3] = Error[i,3]/20
	Error[i,4] = Error[i,4]/20;
	Error[i,5] = Error[i,5]/20;

	#println("Error: ",Error[i,2]);
        i=i+1;
        println((N/MaxIT)*100);
    end
    return Error;
end
function BanchDispersion(r)
	Omega = zeros(100,4);
	for N in 1:1:100
		D=[0.0,1.0];
		h = (D[2]-D[1])/(N+1);
		dt = r*h;
		d1=(-1)*ones(1,N-1);
		d2=[];
		for i in 1:length(d1)+1
			d2=push!(d2,2);
		end
		d2=convert(Array{Float64}, d2);
		A=spdiagm(-1=>d1[1,:],0=>d2,1=>d1[1,:]);
		A=(1/(h^2))*A;

		Ω1 = [];
		Ω2 = [];
		Ω3 = [];	

		for eig in eigvals(Matrix(A))
			ω1 = (2/dt)*asin(sqrt((0.25*dt^2)*eig));
			Ω1 = push!(Ω1,ω1);

			ω2 = sqrt(eig);
			Ω2 = push!(Ω2,ω2);

			ω3 = (1/dt)*acos((1-0.25*eig*dt^2)/(0.25*(eig*dt^2)+1));
			Ω3 = push!(Ω3,ω3);
		end
		#println("Numerical anglar speed: ",Ω[1]);
		Omega[N,1]=h;
		Omega[N,2]=abs(π-Ω1[1]);
		Omega[N,3]=abs(π-Ω2[1]);
		Omega[N,4]=abs(π-Ω3[1]);
	end
	return Omega;
end
function NewmarkBanchTimeWindow(r,T)
	D=[0.0,1.0];
	Ini=[0.0,0.0];
	f(x)=sin(π*x);
	g(x)=0*x;
	e(x,t)=sin(π*x)*cos(π*t)
	MaxIT = 1010;
	Error = zeros(Int64(MaxIT/10),5);
	Energy = [];
    i=1;
    for N in 10:10:MaxIT
        for j in 1:20
            τ = @elapsed begin
                Y,h = NewmarkWave(f,g,N,D,Ini,r,T);
            end
            E = []
	    Amp = []
            x=D[1]:h:D[2];
            t=0;
	    for k in 2:length(Y)-1
		y = Y[k];
                t=t+r*h
		A = maximum(e.(x,t));	
		Amp = push!(Amp,abs(abs(A)-abs(maximum(y))));
                E = push!(E,norm(y-e.(x,t),Inf))
		Energy = push!(Energy,abs(0.5*(π^2)-DirichletEnergy(Y[k-1],Y[k+1],h,r*h))); 
            end
    	    x=D[1]:h:D[2];
            Error[i,1] = Error[i,1]+h;
            Error[i,2] = Error[i,2]+norm(E,Inf);
            Error[i,3] = Error[i,3]+τ;
	    Error[i,4] = Error[i,4]+Amp[end];
	   
	    Error[i,5] = Error[i,5]+Energy[end];
        end
        Error[i,1]=Error[i,1]/20;
        Error[i,2]=Error[i,2]/20;
        Error[i,3] = Error[i,3]/20;
	Error[i,4]=Error[i,4]/20;
	Error[i,5]=Error[i,5]/20;
        i=i+1;
        println((N/MaxIT)*100);
    end
    return Error;
end
function MillerBanchTimeWindow(r,T)
	D=[0.0,1.0];
	Ini=[0.0,0.0];
	f(x)=sin(π*x);
	g(x)=0*x;
	e(x,t)=sin(π*x)*cos(π*t)
    Error = zeros(100,5);
    i=1;
    for N in 1:10:1000
        for j in 1:20
            τ = @elapsed begin
                Y,h = MillerWave(f,g,N,D,Ini,r,T);
            end
            E = []
	    Amp = []
            x=D[1]:h:D[2];
            t=r*h;
            for y in Y
		t=t+(r*h);
		A = maximum(e.(x,t));	
		Amp = push!(Amp,abs(abs(A)-abs(maximum(y))));
                E = push!(E,norm(y-e.(x,t),Inf))
            end
    		x=D[1]:h:D[2];
            Error[i,1] = Error[i,1]+h;
            Error[i,2] = Error[i,2]+norm(E,Inf);
            Error[i,3] = Error[i,3]+τ;
	    Error[i,4] = Error[i,4]+Amp[end];
        end
        Error[i,1]=Error[i,1]/20;
        Error[i,2]=Error[i,2]/20;
        Error[i,3] = Error[i,3]/20;
	Error[i,4]=Error[i,4]/20;
        i=i+1;
        println((N/1000)*100)
    end
    return Error;
end
function MillerBanchTimeStep(T)
	D=[0.0,1.0];
	Ini=[0.0,0.0];
	f(x)=sin(π*x);
	g(x)=0*x;
	e(x,t)=sin(π*x)*cos(π*t)
    Error = zeros(100,3);
    i=1;
	N=10000-1;
	Δx=(D[2]-D[1])/(N+1);
    for r in 1:10:1000
		println((r/1000)*100)
        for j in 1:1
            τ = @elapsed begin
                Y,h = MillerWave(f,g,N,D,Ini,r,T);
            end
            E = []
            x=D[1]:h:D[2];
            t=r*h;
            for y in Y
                t=t+(r*h);
                E = push!(E,norm(y-e.(x,t),Inf))
            end
            Error[i,1] = Error[i,1]+r;
            Error[i,2] = Error[i,2]+r*h;
            Error[i,3] = Error[i,3]+norm(E,Inf);
        end
        Error[i,1]=Error[i,1]/1;
        Error[i,2]=Error[i,2]/1;
        Error[i,3]=Error[i,3]/1;
        i=i+1;
    end
    return Error, Δx;
end
function MillerBanchSpaceStep(T)
	D=[0.0,1.0];
	Ini=[0.0,0.0];
	f(x)=sin(π*x);
	g(x)=0*x;
	e(x,t)=sin(π*x)*cos(π*t)
    Error = zeros(100,3);
    i=1;
	Δt = 10^(-4);
    for N in 1:10:1000
		println((N/1000)*100)
		Δx=(D[2]-D[1])/(N+1);
		r = (Δt/Δx);
        for j in 1:1
            τ = @elapsed begin
                Y,h = MillerWave(f,g,N,D,Ini,r,T);
            end
            E = []
            x=D[1]:h:D[2];
            t=r*h;
            for y in Y
                t=t+(r*h);
                E = push!(E,norm(y-e.(x,t),Inf))
            end
            Error[i,1] = Error[i,1]+N;
            Error[i,2] = Error[i,2]+h;
            Error[i,3] = Error[i,3]+norm(E,Inf);
        end
        Error[i,1]=Error[i,1]/1;
        Error[i,2]=Error[i,2]/1;
        Error[i,3]=Error[i,3]/1;
        i=i+1;
    end
    return Error, Δt;
end
function NewmarkBanchSpaceStep(T)
	D=[0.0,1.0];
	Ini=[0.0,0.0];
	f(x)=sin(π*x);
	g(x)=0*x;
	e(x,t)=sin(π*x)*cos(π*t)
    Error = zeros(100,3);
    i=1;
	N=2000;
	Δt = 10^(-4);
    for N in 1:10:1000
		println((N/1000)*100)
		Δx=(D[2]-D[1])/(N+1);
		r = (Δt/Δx);
        for j in 1:1
            τ = @elapsed begin
                Y,h = NewmarkWave(f,g,N,D,Ini,r,T);
            end
            E = []
            x=D[1]:h:D[2];
            t=0;
            for y in Y
                t=t+r*h
                E = push!(E,norm(y-e.(x,t),Inf))
            end
    		x=D[1]:h:D[2];
            Error[i,1] = Error[i,1]+N;
            Error[i,2] = Error[i,2]+h;
            Error[i,3] = Error[i,3]+norm(E,Inf);
        end
        Error[i,1]=Error[i,1]/1;
        Error[i,2]=Error[i,2]/1;
        Error[i,3]=Error[i,3]/1;
        i=i+1;
    end
    return Error, Δt;
end
function NewmarkBanchTimeStep(T)
	D=[0.0,1.0];
	Ini=[0.0,0.0];
	f(x)=sin(π*x);
	g(x)=0*x;
	e(x,t)=sin(π*x)*cos(π*t)
    Error = zeros(100,3);
    i=1;
	N=10000-1;
	Δx=(D[2]-D[1])/(N+1);
    for r in 1:10:1000
		println((r/1000)*100)
        for j in 1:1
            τ = @elapsed begin
                Y,h = NewmarkWave(f,g,N,D,Ini,r,T);
            end
            E = []
            x=D[1]:h:D[2];
            t=0;
            for y in Y
                t=t+r*h
                E = push!(E,norm(y-e.(x,t),Inf))
            end
    		x=D[1]:h:D[2];
            Error[i,1] = Error[i,1]+r;
            Error[i,2] = Error[i,2]+r*h;
            Error[i,3] = Error[i,3]+norm(E,Inf);
        end
        Error[i,1]=Error[i,1]/1;
        Error[i,2]=Error[i,2]/1;
        Error[i,3]=Error[i,3]/1;
        i=i+1;
    end
    return Error, Δx;
end
function MillerZBanchTimeWindow(r,T)
	D=[0.0,1.0];
	Ini=[0.0,0.0];
	f(x)=sin(π*x);
	g(x)=0*x;
	e(x,t)=sin(π*x)*cos(π*t)
    Error = zeros(100,3);
    i=1;
    for N in 1:10:1000
        for j in 1:20
            τ = @elapsed begin
                Y,h = MillerZWave(f,g,N,D,Ini,r,T);
            end
            E = []
            x=D[1]:h:D[2];
            t=r*h;
            for y in Y
                t=t+(r*h);
                E = push!(E,norm(y-e.(x,t),Inf))
            end
    		x=D[1]:h:D[2];
            Error[i,1] = Error[i,1]+h;
            Error[i,2] = Error[i,2]+norm(E,Inf);
            Error[i,3] = Error[i,3]+τ;
        end
        Error[i,1]=Error[i,1]/20;
        Error[i,2]=Error[i,2]/20;
        Error[i,3]=Error[i,3]/20;
        i=i+1;
        println((N/1000)*100)
    end
    return Error;
end
function SpeedLeapfrogBanchTimeWindow(r,T,v)
	D=[0.0,1.0];
	Ini=[0.0,0.0];
	f(x)=sin(π*x);
	g(x)=0*x;
	e(x,t)=sin(π*x)*cos(π*t*v)
    Error = zeros(100,3);
    i=1;
    for N in 1:10:1000
        for j in 1:20
            τ = @elapsed begin
                Y,h = LeapfrogWave(f,g,N,D,Ini,r,T,v);
            end
            E = []
            x=D[1]:h:D[2];
            t=r*h;
            for y in Y
                t=t+r*h
                E = push!(E,norm(y-e.(x,t),Inf))
            end
    		x=D[1]:h:D[2];
            Error[i,1] = Error[i,1]+h;
            Error[i,2] = Error[i,2]+norm(E,Inf);
            Error[i,3] = Error[i,3]+τ;
        end
        Error[i,1] = Error[i,1]/20;
        Error[i,2] = Error[i,2]/20;
        Error[i,3] = Error[i,3]/20
        i=i+1;
        println((N/1000)*100);
    end
    return Error;
end
function SpeedNewmarkBanchTimeWindow(r,T,v)
	D=[0.0,1.0];
	Ini=[0.0,0.0];
	f(x)=sin(π*x);
	g(x)=0*x;
	e(x,t)=sin(π*x)*cos(π*t*v)
    Error = zeros(100,3);
    i=1;
    for N in 1:10:1000
        for j in 1:20
            τ = @elapsed begin
                Y,h = NewmarkWave(f,g,N,D,Ini,r,T,v);
            end
            E = []
            x=D[1]:h:D[2];
            t=0;
            for y in Y
                t=t+r*h
                E = push!(E,norm(y-e.(x,t),Inf))
            end
    		x=D[1]:h:D[2];
            Error[i,1] = Error[i,1]+h;
            Error[i,2] = Error[i,2]+norm(E,Inf);
            Error[i,3] = Error[i,3]+τ;
        end
        Error[i,1] = Error[i,1]/20;
        Error[i,2] = Error[i,2]/20;
        Error[i,3] = Error[i,3]/20
        i=i+1;
        println((N/1000)*100);
    end
    return Error;
end
function SpeedMillerBanchTimeWindow(r,T,v)
	D=[0.0,1.0];
	Ini=[0.0,0.0];
	f(x)=sin(π*x);
	g(x)=0*x;
	e(x,t)=sin(π*x)*cos(π*t*v)
    Error = zeros(100,3);
    i=1;
    for N in 1:10:1000
        for j in 1:20
            τ = @elapsed begin
                Y,h = MillerWave(f,g,N,D,Ini,r,T,v);
            end
            E = []
            x=D[1]:h:D[2];
            t=r*h;
            for y in Y
                t=t+(r*h);
				E = push!(E,norm(y-e.(x,t),Inf))
            end
    		x=D[1]:h:D[2];
            Error[i,1] = Error[i,1]+h;
            Error[i,2] = Error[i,2]+norm(E,Inf);
            Error[i,3] = Error[i,3]+τ;
        end
        Error[i,1] = Error[i,1]/20;
        Error[i,2] = Error[i,2]/20;
        Error[i,3] = Error[i,3]/20
        i=i+1;
        println((N/1000)*100);
    end
    return Error;
end
function WaveBanch(opt)
	##########
	#[1] Benchmark the Leapfrog method v=1. Showing CFL
	#[2] Benchmark the Miller method v=1;
	#[3] Benchmark the Newmark mathod v=1;
	#[4] Comparison of the Performance for v=1
	#[5] Benchmark and Comparison of the Performance for v=2.8
    r=[0.5,0.75,1.0,1.05,1.5]
    rs=[0.5,0.75,0.99,1.05,1.5]
    if opt==1
        ############| LEAPFROG |##########
        figure()
        for i in 1:5
            R3 = LeapfrogBanchTimeWindow(rs[i],[0.0,π/2]);
            loglog(R3[:,1],R3[:,2],marker="o",label=string("Leapfrog r=",round(r[i];digits=3)))
        end
        title(L"Leapfrog Method Error Evaluated in t$\in [0,\frac{\pi}{2}]$, $c=1$, $||\cdot||_2$")
	ylabel(L"Error $||\cdot||_\infty$")
	xlabel(L"Step Size $h_x$")
        legend(loc=0,borderaxespad=0);
		ylim(10^(-8),10^(10));
		println(R3);
	savefig("A1.png");
        figure()
        for i in 1:3
             R4 = LeapfrogBanchTimeWindow(rs[i],[0,π/2]);
             loglog(R4[:,3],R4[:,2],marker="o",label=string("Leapfrog r=",round(r[i];digits=3)))
        end
        title(L"Leapfrog Method Performance Evaluated in t$\in [0,\frac{\pi}{2}]$, $c=1$, $||\cdot||_2$")
	ylabel(L"Error $||\cdot||_\infty$")
	xlabel("Computational Time in s")
        legend(loc=0,borderaxespad=0);
	savefig("B1.png");
    elseif opt==2
        ##########| MILLER-GRIFFITHS|##########
        figure()
        for i in [1,2,3,5]
            R3 = MillerBanchTimeWindow(r[i],[0,π/2]);
            loglog(R3[:,1],R3[:,2],marker="o",label=string("Miller-Griffiths r=",round(r[i];digits=3)))
        end
        title(L"Miller Griffiths Scheme Error Evaluated in t $\in [0,\frac{\pi}{2}]$, $c=1$, $||\cdot||_\infty$")
	ylabel(L"Error $||\cdot||_\infty$")
	xlabel(L"Step Size $h_x$")
        legend(loc=0,borderaxespad=0)
	savefig("A2.png")
        figure()
        for i in [1,2,3,5]
            R3 = MillerBanchTimeWindow(r[i],[0,π/2]);
            loglog(R3[:,2],R3[:,3],marker="o",label=string("Miller-Griffiths r=",round(r[i];digits=3)))
        end
        title(L"Miller Griffiths Scheme Performance Evaluated in t$\in [0,\frac{\pi}{2}]$, $c=1$, $||\cdot||_\infty$")
        legend(loc=0,borderaxespad=0);
	ylabel(L"Error $||\cdot||_\infty$")
	xlabel("Computational Time in s")
	savefig("B2.png")
    elseif opt==3
        figure()
        for i in [1,2,3,5]
            R3 = NewmarkBanchTimeWindow(r[i],[0,π/2]);
            loglog(R3[:,1],R3[:,2],marker="o",label=string("Newmark r=",round(r[i];digits=3)))
        end
        title(L"Newmark Integration Error Evaluated in t$\in [0,\frac{\pi}{2}]$, $c=1$, $||\cdot||_\infty$")
        legend(loc=0,borderaxespad=0);
	ylabel(L"Error $||\cdot||_\infty$")
	xlabel(L"Step Size $h_x$")
	savefig("A3.png");
        figure()
        for i in [1,2,3,5]
            R3 = NewmarkBanchTimeWindow(r[i],[0,π/2]);
            loglog(R3[:,2],R3[:,3],marker="o",label=string("Newmark r=",round(r[i];digits=3)))
        end
        title(L"Newmark Integration Performance Evaluated in t$\in [0,\frac{\pi}{2}]$, $c=1$, $||\cdot||_\infty$")
        legend(loc=0,borderaxespad=0);
	ylabel(L"Error $||\cdot||_\infty$")
	xlabel("Computational Time in s")
	savefig("B3.png");
    elseif opt==4
        ##########| COMPARED |##########
        figure()
        R1 = LeapfrogBanchTimeWindow(0.98,[0,π/2]);
	R6 = LeapfrogBanchTimeWindow(0.5,[0,π/2]);
        R2= NewmarkBanchTimeWindow(1.0,[0,π/2]);
        R3= NewmarkBanchTimeWindow(5.0,[0,π/2]);
        R5= MillerBanchTimeWindow(5.0,[0,π/2]);
        R4 = MillerBanchTimeWindow(1.0,[0,π/2]);
        loglog(R1[:,3],R1[:,2],marker="o",label=string("Leapfrog r=1"))
	loglog(R6[:,3],R6[:,2],marker="o",label=string("Leapfrog r=0.5"))
        loglog(R4[:,3],R4[:,2],marker="o",label=string("Miller-Griffiths r=1"))
        loglog(R2[:,3],R2[:,2],marker="o",label=string("Newmark r=1"))
        loglog(R5[:,3],R5[:,2],marker="o",label=string("Miller-Griffiths r=5"))
        loglog(R3[:,3],R3[:,2],marker="o",label=string("Newmark r=5"))
        legend(loc=2,borderaxespad=0);
        title(L"Performance Comparison $c=1$, $||\cdot||_\infty$");
	ylabel(L"Error $||\cdot||_\infty$")
	xlabel(L"Computational Time in s")
	savefig("A4.png");
	elseif opt == 5
		figure()
		R1 = SpeedLeapfrogBanchTimeWindow(0.35,[0,π/2],2.8);
		R2= SpeedLeapfrogBanchTimeWindow(0.75,[0,π/2],2.8);
		R3= SpeedLeapfrogBanchTimeWindow(1.0,[0,π/2],2.8);
		loglog(R1[2:100,1],R1[2:100,2],marker="o",label=string("Leapfrog r=0.35"))
		loglog(R2[2:100,1],R2[2:100,2],marker="o",label=string("Leapfrog r=0.75"))
		loglog(R3[2:100,1],R3[2:100,2],marker="o",label=string("Leapfrog r=1"))
		title(L"Leapfrog Method Error Evaluated in t$\in [0,\frac{\pi}{2}]$, $c=2.8m/s$")
		legend(loc=0,borderaxespad=0);
		ylabel(L"Error $||\cdot||_\infty$")
		xlabel(L"Step Size $h_x$")
		savefig("A5.png");
		figure()
		R4 = SpeedNewmarkBanchTimeWindow(0.35,[0,π/2],2.8);
		R5= SpeedNewmarkBanchTimeWindow(0.75,[0,π/2],2.8);
		R6= SpeedNewmarkBanchTimeWindow(1.0,[0,π/2],2.8);
		loglog(R4[2:100,1],R4[2:100,2],marker="o",label=string("Newmark r=0.35"))
    	loglog(R5[2:100,1],R5[2:100,2],marker="o",label=string("Newmark r=0.75"))
		loglog(R6[2:100,1],R6[2:100,2],marker="o",label=string("Newmark r=1"))
		title(L"Newmark Method Error Evaluated in t$\in[0,\frac{\pi}{2}]$, $c=2.8m/s$")
		legend(loc=0,borderaxespad=0);
		ylabel(L"Error $||\cdot||_\infty$")
		xlabel(L"Step Size $h_x$")
		savefig("B5.png");
		figure()
		R7 =SpeedMillerBanchTimeWindow(0.35,[0,π/2],2.8);
		R8= SpeedMillerBanchTimeWindow(0.75,[0,π/2],2.8);
		R9= SpeedMillerBanchTimeWindow(1.0,[0,π/2],2.8);
		loglog(R7[2:100,1],R7[2:100,2],marker="o",label=string("Miller-Griffiths r=0.35"))
    	loglog(R8[2:100,1],R8[2:100,2],marker="o",label=string("Miller-Griffiths r=0.75"))
		loglog(R9[2:100,1],R9[2:100,2],marker="o",label=string("Miller-Griffiths r=1"))
		title(L"Miller-Griffiths Method Error Evaluated in t$\in [0,\frac{\pi}{2}]$, $c=2.8m/s$")
		legend(loc=0,borderaxespad=0);
		ylabel(L"Error $||\cdot||_\infty$")
		xlabel(L"Step Size $h_x$")
		savefig("C5.png");
	elseif opt == 6
		figure()
		R1 = SpeedLeapfrogBanchTimeWindow(0.35,[0,π/2],2.8);
		R2 = SpeedMillerBanchTimeWindow(1,[0,π/2],2.8);
		R3 = SpeedMillerBanchTimeWindow(5,[0,π/2],2.8);
		R4 = SpeedNewmarkBanchTimeWindow(0.75,[0,π/2],2.8);
		R5 = SpeedNewmarkBanchTimeWindow(1,[0,π/2],2.8);
		loglog(R1[2:100,2],R1[2:100,3],marker="o",label=string("Leapfrog r=0.35"))
    	loglog(R2[2:100,2],R2[2:100,3],marker="o",label=string("Miller-Griffiths r=1"))
		loglog(R3[2:100,2],R3[2:100,3],marker="o",label=string("Miller-Griffiths r=5"))
    	loglog(R4[2:100,2],R4[2:100,3],marker="o",label=string("Newmark r=0.75"))
		loglog(R5[2:100,2],R5[2:100,3],marker="o",label=string("Newmark r=1"))
		title(L"Performance Comparison in t$\in[0,\frac{\pi}{2}]$, $c=2.8m/s$, $||\cdot||_\infty$")
		ylabel(L"Error $||\cdot||_\infty$")
		xlabel("Computational Time in s")
                legend(loc=0,borderaxespad=0);
		savefig("A6.png");
		legend(loc=0,borderaxespad=0);

	elseif opt == 7
		figure()
		R1 = SpeedLeapfrogBanchTimeWindow(0.1,[0,π/2],5);
		R2 = SpeedMillerBanchTimeWindow(0.5,[0,π/2],5);
		R3 = SpeedMillerBanchTimeWindow(1,[0,π/2],5);
		R4 = SpeedNewmarkBanchTimeWindow(0.5,[0,π/2],5);
		R5 = SpeedNewmarkBanchTimeWindow(1,[0,π/2],5);
		loglog(R1[2:100,2],R1[2:100,3],marker="o",label=string("Leapfrog r=0.1"))
    	loglog(R2[2:100,2],R2[2:100,3],marker="o",label=string("Miller-Griffiths r=0.5"))
		loglog(R3[2:100,2],R3[2:100,3],marker="o",label=string("Miller-Griffiths r=1"))
    	loglog(R4[2:100,2],R4[2:100,3],marker="o",label=string("Newmark r=0.5"))
		loglog(R5[2:100,2],R5[2:100,3],marker="o",label=string("Newmark r=1"))
		title(L"Performance Comparison in t$\in[0,\frac{\pi}{2}]$, $c=5m/s$, $||\cdot||_\infty$")
		legend(loc=0,borderaxespad=0);
		ylabel(L"Error $||\cdot||_\infty$")
		xlabel("Computational Time in s")
		savefig("A7.png");
	elseif opt==8
		##########| MILLER-GRIFFITHS MODIFICATO|##########
		figure()
		for i in 1:5
			R3 = MillerZBanchTimeWindow(r[i],[0,π/2]);
			loglog(R3[:,1],R3[:,2],marker="o",label=string("Miller-Griffiths r=",round(r[i];digits=1)))
		end
		title(L"Miller Griffiths Scheme Error Evaluated in t=$ [0,\frac{\pi}{2}]$, $v=1$, $||\cdot||_\infty$")
		legend(loc=0,borderaxespad=0);
		figure()
		for i in 1:5
			R3 = MillerZBanchTimeWindow(r[i],[0,π/2]);
			loglog(R3[:,2],R3[:,3],marker="o",label=string("Miller-Griffiths r=",round(r[i];digits=1)))
		end
		title(L"Miller Griffiths Scheme Performance Evaluated in t=$ [0,\frac{\pi}{2}]$, $v=1$, $||\cdot||_\infty$")
		legend(loc=0,borderaxespad=0);
	elseif opt==9
				
		figure()
		R1, Δx = MillerBanchTimeStep([0.0,π/2]);
		loglog(R1[:,2],R1[:,3],marker="o",label=L"$h_x=0.0001$");
		R2, Δx = MillerBanchSpaceStep([0.0,π/2]);
		loglog(R2[:,2],R2[:,3],marker="o",label=L"$h_t=0.0001$");
		ylabel(L"Error $ ||\cdot||_\infty $");
		title("Miller-Griffiths Scheme");
		legend(loc=0,borderaxespad=0);
		savefig("A9.png");
		return R1, R2;
	elseif opt==10
		figure()
		figure()
		R1, Δx = NewmarkBanchTimeStep([0.0,π/2]);
		loglog(R1[:,2],R1[:,3],marker="o",label=L"$h_x=0.0001$");
		ylabel(L"Error $ ||\cdot||_\infty $");
		title("Newmark Integration");

		R2, Δx = NewmarkBanchSpaceStep([0.0,π/2]);
		loglog(R2[:,2],R2[:,3],marker="o",label=L"$h_t=0.0001$");
		legend(loc=0,borderaxespad=0);
		savefig("A10.png");
		return R1, R2;
	elseif opt==11
        	for i in 1:3
			figure()
			R3 = LeapfrogBanchTimeWindow(r[i],[0.0,π/2]);
			R1 = MillerBanchTimeWindow(r[i],[0.0,π/2]);
			R2 = NewmarkBanchTimeWindow(r[i],[0.0,π/2]);
            		loglog(R3[:,1],R3[:,4],marker="o",label=string("Leapfrog r=",round(r[i];digits=1)))
        	
            		loglog(R2[:,1],R2[:,4],marker="o",label=string("Newmark r=",round(r[i];digits=1)))
											    
            		loglog(R1[:,1],R1[:,4],marker="o",label=string("Miller-Griffiths r=",round(r[i];digits=1)))
											    
        		title(string(L"Dissipation Evaluated in t$\in  [0,\frac{\pi}{2}]$, $c=1$,$r=",round(rs[i];digits=3),L"$, $||\cdot||_2$"));
        		legend(loc=0,borderaxespad=0);
			savefig("A11.png");
		end
	elseif opt==12
		for i in 1:3
			figure()
			R1 = BanchDispersion(r[i]);
			loglog(R1[:,1],R1[:,2],marker="o",label=string("Leapfrog r=",round(r[i];digits=1)));
			loglog(R1[:,1],R1[:,3],marker="o",label=string("Newmark r=",round(r[i];digits=1)));
			loglog(R1[:,1],R1[:,4],marker="o",label=string("Miller-Griffiths r=",round(r[i];digits=1)));
			legend(loc=0,borderaxespad=0);
			title(string(L"Dispersion, $r=",round(rs[i];digits=3),L"$"));
			xlabel(L"Spatial Step Size $h_x$")
			ylabel(L"Angular Speed Error $e_\omega$");
			legend(loc=0,borderaxespad=0);
			savefig("A12.png");
		end
      elseif opt==13	
        figure()
        for i in 1:3
             R5 = LeapfrogBanchTimeWindow(r[i],[0,π/2]);
	     saveArray(string("Leapfrog_",r[i],".csv"),["h","Error","Timing", "Amplitude","Energy"],R5);
             loglog(R5[:,1],R5[:,5],marker="o",label=string("Leapfrog r=",round(r[i];digits=1)))
        end
        title(L"Leapfrog Method Energy Evaluated in t$\in [0,\frac{\pi}{2}]$, $c=1$, $||\cdot||_2$")
        legend(loc=0,borderaxespad=0);
	xlabel(L"Space Step Size $h_x$");
	ylabel(L"Energy Error $e_E$");
	savefig("A13.png");
	return R5;
      elseif opt==14
        figure()
        for i in 1:3
             R6 = NewmarkBanchTimeWindow(r[i],[0,π/2]);
	     saveArray(string("Newmark_",r[i],".csv"),["h","Error","Timing", "Amplitude","Energy"],R6);
             loglog(R6[:,1],R6[:,5],marker="o",label=string("Newmark r=",round(r[i];digits=1)))
        end
        title(L"Newmark Method Energy Evaluated in t$\in [0,\frac{\pi}{2}]$, $c=1$, $||\cdot||_2$")
        legend(loc=0,borderaxespad=0);
	legend(loc=0,borderaxespad=0);
	xlabel(L"Space Step Size $h_x$");
	ylabel(L"Energy Error $e_E$");
	savefig("A14.png");
    	return R6;
      elseif opt==15
        figure()
        for i in 1:3
             R6 = MillerBanchTimeWindow(r[i],[0,π/2]);
	     saveArray(string("Miller_",r[i],".csv"),["h","Error","Timing", "Amplitude","Energy"],R6);
             loglog(R6[:,1],R6[:,5],marker="o",label=string("Newmark r=",round(r[i];digits=1)))
        end
        title(L"Miller Method Energy Evaluated in t$\in [0,\frac{\pi}{2}]$, $c=1$, $||\cdot||_2$")
        legend(loc=0,borderaxespad=0);
	legend(loc=0,borderaxespad=0);
	xlabel(L"Space Step Size $h_x$");
	ylabel(L"Energy Error $e_E$");
	savefig("A15.png");
    	return R6;
      end
end
