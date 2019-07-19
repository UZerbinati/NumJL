using PyPlot
function FiniteDirichletExample()
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
