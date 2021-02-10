function matrix_construct(N::Int, qf::Float64)
    q = collect(0:0.1:qf)*1im
    d = [4r^2 for r in 0:N]
    l = ones(N-1)
    l[1] = sqrt(2)
    A = spdiagm(0 => d, 1 => q[end].*l, -1=> q[end].*l)
    return A
end
