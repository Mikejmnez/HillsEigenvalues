function matrix_construct(N::Int64, q::Complex{Float64}, alphas::Array{Float64, 1})
    d = [4r^2 for r in 0:N] # N+1 array
    A = spdiagm(0 => d)
    l = ones(N)
    if length(alphas) == 2  # This is Mathieu's canonical form
        l[1] = sqrt(2)
    end
    for i in 1:length(alphas[2:end])+1
        as = alphas[i]*ones(N-i+1)
        A = A + spdiagm(i => q.*as, -i=> q.*as)
    end
    return A
end

function eigval_iterator(N::Int64, q::Array{Complex{Float64}, 1}, alphas::Array{Float64, 1})
    if q[1] == 0
        q[1] = 0.5*(q[1]+q[2]);
    end
    A = matrix_construct(N, q[1], alphas);
    l, = eigs(A, nev= N-1, which=:SM);
    l = round.(l, digits=2);
    I = sortperm(l, lt = (x,y) -> real(x)==real(y) ? imag(x)<imag(y) : real(x)<real(y));
    a = l[I];
    for qf in q[2:end]
        A = matrix_construct(N, qf, alphas);
        l, = eigs(A, nev=N-1, which=:SM);
        l = round.(l, digits=2);
        I = sortperm(l, lt = (x,y) -> real(x)==real(y) ? imag(x)<imag(y) : real(x)<real(y));
        a = hcat(a, l[I]);
    end
    vals = Dict("a0" => a[1, :], "q"=> q);
    for n in 1:length(l)-1
        vals["a"*string(2n)] = a[n+1, :];
    end
    return vals
end
