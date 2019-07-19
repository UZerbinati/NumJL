function bisection(f,x,tol,itmax=10^6)
	i=0;
	while(i<itmax)
		c=(x[1]+x[2])/2;
		if(abs(f(c))<tol)
			return c;
		end
		if(f(c)*f(x[2])<0)
			x[1]=c;
		else
			x[2]=c;
		end
		i=i+1;
	end
end
function newton(f,x0,tol=0.0001,itmax=10^6)
#newton, is a function that applis the newton method near by a point x0 to find closest zero, of a function f(x)=x^2-4
#Example: newton(f,5);
	x=[Float64(x0)];
	for i in 1:1:itmax
		y=x[i]-(f(x[i])/diff(f,[x[i]])[1]);
		x=append!(x,y);
		if abs(f(y))<tol
			return y;
		end
	end
end
function complexNewton(f,df,x0,tol=0.0000001,itmax=10^5)
	x=[x0];
	y=[];
	for i in 1:1:itmax
		y=x[i]-(f(x[i])/df(x[i]));
		x=append!(x,y);
		if abs(f(y))<tol
			return y;
		end
	end
	return y;
end
function Newton2D(F,DF,x0,tol=0.00000001,itmax=10^6)
	x=[];
	x=push!(x,x0);
	y=[];
	for i in 1:1:itmax
		#println(x[i])
		#println(inv(DF(x[i][1],x[i][2]))*F(x[i][1],x[i][2]))
		y=x[i]-inv(DF(x[i][1],x[i][2]))*F(x[i][1],x[i][2]);
		#println(y)
		#println(sum(abs.(F(y[1],y[2]))))
		x=push!(x,y);
		if sum(abs.(F(y[1],y[2])))<tol
			return y;
		end
	end
	return y;
end
