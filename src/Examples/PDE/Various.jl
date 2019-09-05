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
