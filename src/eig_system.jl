function matrix_construct(N::Int64, q::Complex{Float64})
    d = [4r^2 for r in 0:N]
    l = ones(N)
    l[1] = sqrt(2)
    A = spdiagm(0 => d, 1 => q.*l, -1=> q.*l)
    return A
end


function eigval_iterator(N::Int64, q::Array{Complex{Float64}, 1})
    if q[1] == 0
        q[1] = 0.5*(q[1]+q[2]);
    end
    A = matrix_construct(N, q[1]);
    l, = eigs(A, nev= N-1, which=:SM);
    l = round.(l, digits=2);
    I = sortperm(l, lt = (x,y) -> real(x)==real(y) ? imag(x)<imag(y) : real(x)<real(y));
    a = l[I];
    for qf in q[2:end]
        A = matrix_construct(N, qf);
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
