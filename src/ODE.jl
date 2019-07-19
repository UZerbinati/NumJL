################################|EULER METHOD|##################################
function eulerMethod(du,x,h,u0)
#eulerMethod, is a function that applies an euler method to a differential equation du, for example du(x)=x, and return a sequence that approximate the solution of the Cauchy problem:
#du(x)=x and du(0)=u0
#Example: eulerMethod(du,[-5,5],h=0.1,1)
#where h is the step lengh for the Euler method.
	t=minimum(x):h:maximum(x);
	eta=[]
	eta = append!(eta,u0);
	for i in 2:length(t)
		it =  eta[i-1]+h*du(eta[i-1],t[i-1]);
		eta=append!(eta,it);
	end
	return eta,t;
end
function eulerImplicitMethod(du,x,h,u0)
#eulerImplicitMethod, is a function that applies an implicit euler method to a differential equation du, for example du(x)=x, and return a sequence that approximate the solution of the Cauchy problem:
#du(x)=x and du(0)=u0
#Example: eulerImplicitMethod(du,[-5,5],h=0.1,1)
#where h is the step lengh for the Euler method.
	t=minimum(x):h:maximum(x);
	eta=[u0];
	for i in 2:length(t)
		g(s) =  eta[i-1]+h*du(s,t[i])-s;
		it = newton(g,1);
		eta=append!(eta,it);
	end
	return(eta);
end
function complexEulerImplicitMethod(du,ddu,x,h,u0)
#eulerImplicitMethod, is a function that applies an implicit euler method to a differential equation du, for example du(x)=x, and return a sequence that approximate the solution of the Cauchy problem:
#du(x)=x and du(0)=u0
#Example: eulerImplicitMethod(du,[-5,5],h=0.1,1)
#where h is the step lengh for the Euler method.
	t=minimum(x):h:maximum(x);
	eta=[u0];
	eta2=[abs(u0)];
	for i in 2:length(t)
		g(s) =  eta[i-1]+h*du(s,t[i])-s;
		dg(s) = h*ddu(s,t[i])-(1.0+0.0im);
		it = complexNewton(g,dg,0.1+0.1im);
		eta=append!(eta,it);
		eta2=append!(eta2,abs(it));
	end
	return eta2,t ;
end
function thetaMethod(du,x,theta,h,u0)
#thetaMethod, is a function that applies an implicit the theta method to a differential equation du, for example du(x)=x, and return a sequence that approximate the solution of the Cauchy problem:
#du(x)=x and du(0)=u0
#Example: thetaMethod(du,[-5,5],0.5,h=0.1,1)
#where h is the step lengh for the Euler method.
	t=minimum(x):h:maximum(x);
	eta=[u0];
	for i in 2:length(t)
		g(s) =  eta[i-1]+h*((1-theta))*du(eta[i-1],t[i-1])+h*(theta*du(s,t[i]))-s;
		it = newton(g,1);
		eta=append!(eta,it);
	end
	return(eta);
end
function euler2D(du1,du2,x,h,u0)
#euler2D, is a function that applies the euler method to a system of 2 differential equation.
#Example: "Lotka-Volterra"
#du1(x,y,t)=8*x-2*x*y;
#du2(x,y,t)=-5*y+3*x*y;
#y=euler2D(du1,du2,[0.0,20.0],0.000001,[3.0,4.0])
	t=minimum(x):h:maximum(x);
	y1=[u0[1]];
	y2=[u0[2]];
	for i in 2:length(t)
		it1=y1[i-1]+h*du1(y1[i-1],y2[i-1],t[i-1]);
		it2=y2[i-1]+h*du2(y1[i-1],y2[i-1],t[i-1]);
		y1=push!(y1,it1);
		y2=push!(y2,it2);
	end
	return [y1,y2], t;
end
################################|ADAPTIVE METHOD IVP|###########################
function RK45(du,x,u0,h=0.01,tol=eps())
#This function applies the runghe-kutta-felhberg method for IVP problem.
#Example:
#du(x,t)=x
#y,D=RK45(du,[0.0,5.0],1.0);
#plot(D,y);
	f(t,X)=du(X,t);
	#Compute the constants once
	c30 = 3/8;
	c31 = 3/32;
	c32 = 9/32;
	c40 = 12/13;
	c41 = 1932/2197;
	c42 = -7200/2197;
	c43 = 7296/2197;
	c51 = 439/216;
	c52 = -8;
	c53 = 3680/513;
	c54 = -845/4104;
	c61 = -8/27;
	c62 = 2;
	c63 = -3544/2565;
	c64 = 1859/4104;
	c65 = -11/40;
	cw1 = 25/216;
	cw3 = 1408/2565;
	cw4 = 2197/4104;
	cw5 = -1/5;
	cz1 = 16/135;
	cz3 = 6656/12825;
	cz4 = 28561/56430;
	cz5 = -9/50;
	cz6 = 2/55;
	ce1 = 1/360;
	ce3 = -128/4275;
	ce4 = -2197/75240;
	ce5 = 1/50;
	ce6 = 2/55;
	#Absolute Tollerance
	atol = 1e-13;
	alpha = 0.8;
	k = 0;
	#Initial Time Moment
	i=1;
	tt=[x[1]];
	t=x[1];
	y=[u0];
	wi=u0;
	#If it is the last iteration, then lastit = 1, otherwise lastit = 0
	lastit=0;
	while lastit==0
		#Stretch the step if within 10% of b-t
		if t+1.1*h>x[2]
			h=x[2]-t;
			lastit=1;
		end
		#Compute the step
		s1 = f(t,wi);
		s2 = f(t + 0.25 * h, wi + 0.25*h*s1);
		s3 = f(t + c30 * h, wi + c31 * h * s1 + c32 * h * s2);
		s4 = f(t + c40 * h, wi + c41 * h * s1 + c42 * h * s2 + c43 * h * s3);
		s5 = f(t + h, wi + c51 * h * s1 + c52 * h * s2 + c53 * h * s3 + c54 * h * s4);
		s6 = f(t + 0.5 * h, wi + c61 * h * s1 + c62 * h * s2 + c63 * h * s3 + c64 * h * s4 + c65 * h * s5);
	    	w = wi + h * (cw1 * s1 + cw3 * s3 + cw4 * s4 + cw5 * s5);
	    	z = wi + h * (cz1 * s1 + cz3 * s3 + cz4 * s4 + cz5 * s5 + cz6 * s6);
	    	e = h * norm(ce1 * s1 + ce3 * s3 + ce4 * s4 + ce5 * s5 + ce6 * s6);
		#Target Tolerance For this Step
		T=tol*norm(wi)+atol;
		if e<= T
			t=t+h;
			h=alpha*h*(T/e)^0.2;
			i=i+1;
			tt=append!(tt,t);
			wi=z;
			y=append!(y,z);
			k=0;
		elseif k==0
			h=alpha*h*(T/e)^0.2;
			k=k+1;
			lastit=0;
		else
			h=h/2;
			lastit=0;
		end
	end
	return y, tt;
end
############################################|BVP|###############################
function shootingBVP(du1,du2,D,u0,init,tol=eps(),h=0.000001)
#ShootingBVP, is a function that solve a second derivative porblem  with boundery condition. You need to convert the second order differential equation into a system of differential equation.
#Example:
#du1(x,y,t)=y;
#du2(x,y,t)=9.8;
#y=shootingBVP(du1,du2,[0.0,5.0],[0.0,40.0],[1.0,50.0],0.001)
	opt="bisection";
	if (opt=="bisection")
		err=1;
		y3=[];
		alpha=u0[1];
		beta=u0[2];
		#Computing solution differential equation
		y1, t = euler2D(du1,du2,D,h,[alpha,init[1]]);
		y2, t = euler2D(du1,du2,D,h,[alpha,init[2]]);
		y1[1]=y1[1]-beta;
		y2[1]=y2[1]-beta;
		shoot1=init[1];
		shoot2=init[2];
		while(abs(err)>tol)

			c=(shoot1+shoot2)/2;
			#print(c);
			y3,t = euler2D(du1,du2,D,h,[alpha,c]);
			err=y3[1][length(y3[1])]-beta;
			if (y1[1][length(y1[1])]*err)<0
				shoot2=c;
			else
				shoot1=c;
			end

		end
	end
	return y3, t;

end
################################################################################
function leapfrogIntegration(F,x,h,u0)
	#Leapfrog integration solve a problem such as y''=F(y)
	t = x[1]:(h):x[2];
	v = [u0[2]]
	x = [u0[1]]
	for i in 2:length(t)
		v = push!(v, v[i-1]+h*F(x[i-1]));
		x = push!(x, x[i-1]+h*v[i]);
	end
	return x, t;
end
