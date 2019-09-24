function driftExample()
    f0(x) = if (-1 ≤ x ≤  1) return 0.5 else return 0 end
    D = [-4,4];
    x = D[1]:0.01:D[2];
    T = 0.0;
    while(true)
        PyPlot.clf();
        title(string(L"Drift Equation Evolution $\mu =$",mean(drift.(f0,x,T))))
        plot(x,drift.(f0,x,T));
        T = T+0.1;
        msg=readline();
        if msg=="close"
            return;
        end
    end
end
function transportExample()
    f0(x) = if (-1 ≤ x ≤  1) return 0.5 else return 0 end
    D = [-4,4];
    x = D[1]:0.01:D[2];
    T = 0.0;
    while(true)
        PyPlot.clf();
        plot(x,transport.(f0,1,x,T));
        T = T+0.1;
        msg=readline();
        if msg=="close"
            return;
        end
    end
end
function FT()
	N = 2^14 - 1 
	# Sample period
	Ts = 1 / (1.1 * N) 
	# Start time 
	t0 = 0 
	tmax = t0 + N * Ts
	# time coordinate
	t = t0:Ts:tmax

	# signal 
	signal = sin.(2π * 60 .* t) # sin (2π f t) 

	#Fourier Transform of it 
	F = fft(signal) |> fftshift
	freqs = fftfreq(length(t), 1.0/Ts) |> fftshift

	# plots 
	plot(t, signal);
	title("Signal");
	plot(freqs, abs.(F));
	title("Spectrum");
end
function DirichletEnergy(U,U2,h,δ)
	X = 0:h:1;
	
	#U = u(X,0);
	#U2 = u(X,2*δ);
	d1 = (-1)*ones(length(X)-3);
	d2 = ones(length(X)-3);
	#println(d1);
	#println(length(X));
	#println(h);
	D = spdiagm(-1 => d1, 1 => d2);
	∂u = [];
	for i in 1:length(X)
		∂u = push!(∂u, U2[i]-U[i]);
	end
	∂u = (1/(2*δ))*∂u;
	#println(size(U),size(D));
	dx = ((1/(2*h))*D*U[2:end-1]).^2;
	I = TrapezInt(dx,h);
	I = I+TrapezInt((∂u[2:end-1]).^2,δ);
	return I;

end
function Energy()
	u(x,t) = sin.(π*x)*cos.(π*t);
	dx(x,t) = π*cos.(π*x)*cos.(π*t);
	dt(x,t) = -π*sin.(π*x)*sin.(π*t);
	h=0.0001;
	δ = 0.0001;
	X = 0:h:1;
	
	U = u(X,0);
	U1 = u(X,δ);
	U2 = u(X,2*δ);
	
	D = diagm(-1 => (-1)*ones(length(X)-3)[1,:], 1 => ones(length(X)-3));
	println(size(U));
	println(size(D));
	#println(D);
	plot(X[2:end-1],(1/(2*h))*D*U[2:end-1],"r*");
	plot(X,dx(X,0));
	∂u = [];
	for i in 1:length(X)
		∂u = push!(∂u, U2[i]-U[i]);
	end
	∂u = (1/(2*δ))*∂u;
	plot(X,dt(X,δ));
	plot(X,∂u,"r*");
	#I=TrapezInt(U[2:end],h);
	I = TrapezInt(((1/(2*h))*D*U[2:end-1]).^2,h);
	I = I+TrapezInt((∂u[2:end-1]).^2,δ);
	println("Numerical Energy: ",I);
	println("Exact Energy: ", 0.5*(π^2));
end
