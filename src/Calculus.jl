function diff(f,x)
#Recives as input a function f, such as f(x,y)=x^2-4, and the linspace over wich we perform the differentiation.
#Example: diff(f,0:0.1:3)
	x=convert(Array{Float64}, x);
	y=[];
	for i in 1:length(x)
		if x[i] == 0
			x[i]=eps();
		end
		h=sqrt(eps())*x[i];
		xph=x[i]+h;
		dx=xph-x[i];
		y=append!(y,(f(xph)-f(x[i]))/dx);
	end
	return y;
end
function drift(f0,x,t)
	return exp(t)*f0(exp(t)*x);
end
function transport(f0,c,x,t);
	return f0(x-c*t);
end
