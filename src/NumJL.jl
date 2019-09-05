module NumJL
	using LinearAlgebra
	using SparseArrays
	using Dierckx
	
	#For Dumping Testing Only
	using FFTW
	using DSP

	#Only used in some specific case, methods should work even without.
	using PyPlot

	#BASE
	include("LA.jl")
	include("Calculus.jl");
	include("Zeros.jl");
	include("ODE.jl");
	include("PDE.jl");
	include("Geometry.jl");

	#GPU ACCELERATIOn

	#EXAMPLE
	include("Examples/ODE.jl")
	include("Examples/banch.jl")
	include("Examples/PDE/Wave.jl")
	include("Examples/PDE/Heat.jl")
	include("Examples/PDE/FE.jl")
	include("Examples/PDE/FD.jl")
	include("Examples/PDE/Various.jl")
end
