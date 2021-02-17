function matrix_construct(N::Int64, q::Complex{Float64}, alphas::Array{Float64, 1}; flag::String="False", bc::String = "Neumann")
    l = ones(N);
    if bc == "Neumann"
        d = [4r^2 for r in 0:N]; # N+1 array
        if flag == "True"  # This is Mathieu's canonical form. Or if first element is non-zero
            l[1] = sqrt(2) * l[1];
        end
    elseif bc == "Dirichlet"
        d = [4r^2 for r in 1:N+1]; # N+1 array
    else
        error("BC not recognized. Try `Neumann` or `Dirichlet`")
    end
    A = spdiagm(0 => d);
    for i in 0:length(alphas[2:end])
        if alphas[i+1] != 0
            as = alphas[i+1]*l[1:N-i];
            A = A + spdiagm((i+1) => q.*as, -(i+1)=> q.*as);
        end
    end
    return A
end

function eigval_iterator(N::Int64, q::Array{Complex{Float64}, 1}, alphas::Array{Float64, 1}; flag::String="False", bc::String="Neumann")
    if q[1] == 0
        q[1] = 0.5*(q[1]+q[2]);
    end
    A = matrix_construct(N, q[1], alphas; flag, bc);
    l, = eigs(A, nev= N-1, which=:SM);
    l = round.(l, digits=2);
    I = sortperm(l, lt = (x,y) -> real(x)==real(y) ? imag(x)<imag(y) : real(x)<real(y));
    a = l[I];
    for qf in q[2:end]
        A = matrix_construct(N, qf, alphas; flag, bc);
        l, = eigs(A, nev=N-1, which=:SM);
        l = round.(l, digits=2);
        I = sortperm(l, lt = (x,y) -> real(x)==real(y) ? imag(x)<imag(y) : real(x)<real(y));
        a = hcat(a, l[I]);
    end
    if bc == "Neumann"
        var = "a"  # eigenvalue a_2n
        m = 0  # sets lowest index in a
    elseif bc == "Dirichlet"
        var = "b"  # eigenvalue b_2n+2
        m = 2  # sets lowest index in b
    end
    vals = Dict(var*string(m) => a[1, :], "q"=> q);  # initialized Dict
    for n in 1:length(l)-1
        vals[var*string(2n+m)] = a[n+1, :];
    end
    return vals
end
