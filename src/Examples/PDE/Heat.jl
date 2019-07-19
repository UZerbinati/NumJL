function heatExample()
	f(x)=cos(x);
	D=[-3*pi/2,3*pi/2];
	Ini=[0.0,0.0]
	N=50;
	t = 0;
	Δt = 0.5;
	while true
		sleep(0.1)
		PyPlot.clf();
		ax = gca()
		ax[:set_ylim]([-0.3,0.3])
		ax[:set_ylim](-1.1,1.1)
		Y,h = FDTimeStepHeat(f,D,Ini,N,1,t);
		x=D[1]:h:D[2];
		t=t+Δt;
		plot(x,Y)
	end
end
