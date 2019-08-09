######################| Linear Algebra|#################
ndgrid(v::AbstractVector) = copy(v)

function ndgrid(v1::AbstractVector{T}, v2::AbstractVector{T}) where T
    m, n = length(v1), length(v2)
    v1 = reshape(v1, m, 1)
    v2 = reshape(v2, 1, n)
    (repeat(v1, 1, n), repeat(v2, m, 1))
end

function ndgrid_fill(a, v, s, snext)
    for j = 1:length(a)
        a[j] = v[div(rem(j-1, snext), s)+1]
    end
end

function ndgrid(vs::AbstractVector{T}...) where T
    n = length(vs)
    sz = map(length, vs)
    out = ntuple(i->Array{T}(uninitialized, sz), n)
    s = 1
    for i=1:n
        a = out[i]::Array
        v = vs[i]
        snext = s*size(a,i)
        ndgrid_fill(a, v, s, snext)
        s = snext
    end
end

function CGExample()
	#
	# 7 1 0 0 0
	# 1 7 1 0 0
	# 0 1 0.2 1 0
	# 0 0 1 0.1 1
	# 0 0 0 1 0.05
	#
	A = [7 1  0 0 0; 1 7 1 0 0; 0 1 0.3 1 0; 0 0 1 0.2 1; 0 0 0 1 0.1];
	println("----|Matrix A|----");
	println(A);
	println("----|Vector b|----");
	b = transpose([1 2 7 11 3]);

	println("Rank: ",rank(A));
	println("Conditioning number: ", cond(A));
	xe = A\b;
	println("Solution: ",xe);

	#LOW RANK APPROXIMATION
        F = svd(A);
	println("Valori singolari: ",F.S);
	m=4;
	App = F.U[:,1:m] * Diagonal(F.S)[1:m,1:m] * F.Vt[1:m,:];
	println("----| LOW RANK APPROXIMATION OF A|----");
	println(App);	
	
	
	println("----|CG|----");
	x = [];
	x = push!(x,[1.0,1.0,1.0,1.0,1.0]);
	xp = [];
	xp = push!(xp,[1.0,1.0,1.0,1.0,1.0]);

	r = [];
	r = push!(r,b-A*x[end]);
	rp = [];
	rp = push!(rp,b-App*xp[end]);
	
	p = [];
	p = push!(p,r[1]);
	pp = [];
	pp = push!(pp,rp[1]);

	alpha = [];
	alphap = [];
	beta = [];
	betap = [];

	ex = [];
	ep = [];
	ea = [];

	for k in 1:5
		alpha = push!(alpha,(transpose(p[end])*r[end])/(transpose(p[end])* A *p[end]));
		alphap = push!(alphap,(transpose(pp[end])*rp[end])/(transpose(pp[end])* App *pp[end]));
		#println(alpha[end][1]," ",p[end]);
		x = push!(x, x[end] + alpha[end][1]*p[end]);
		xp = push!(xp, xp[end] + alphap[end][1]*pp[end]);

		r = push!(r, b - A*x[end]);
		rp = push!(rp, b - App*xp[end]);
		
		beta = push!(beta,(transpose(p[end]) * A *r[end])/(transpose(p[end])* A *p[end]));
		betap = push!(betap,(transpose(pp[end]) * App *rp[end])/(transpose(pp[end])* App *pp[end]));

		p = push!(p,r[end]-beta[end][1]*p[end]);
		pp = push!(pp,rp[end]-betap[end][1]*pp[end]);
		println("Norm p: ", norm(p[end])," Norm p of low rank approximation: ", norm(pp[end]));

		ex = push!(ex, norm(xe-x[end]));
		ep = push!(ep, norm(p[end]-pp[end]));
		ea = push!(ea, abs(alpha[end][1]-alphap[end][1])/alpha[end][1]);
		
	end
	figure(1)
	plot(ex);
	title("Error of the Solution");
	
	figure(2);
	plot(ep);
	title("Error for Low Rank Approximation Search Direction");
	
	figure(3);
	plot(ea);
	title("Error for Low Rank Approximation Search Step");
	return [[p,pp]]
end
