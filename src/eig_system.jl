function matrix_construct(N::Int, q::Float64)
    d = [4r^2 for r in 0:N]
    l = ones(N)
    l[1] = sqrt(2)
    A = spdiagm(0 => d, 1 => q.*l, -1=> q.*l)
    return A
end
