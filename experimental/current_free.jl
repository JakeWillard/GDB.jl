

function x(theta, epsilon, delta, kappa)

end

function y(theta, epsilon, delta, kappa)

end


function OuterBarrier(N, epsilon, delta, kappa)

    vert = zeros(Float64, (2, N))
    thetas = LinRange(0, -2*pi, N)
    for i=1:N
        vert[1,i] = x(..)
        vert[2,i] = y(...)
    end

    return PolyBarrier(vert)
end


function make_grid(N, ...)

    psi, bx, by = generate_flux(...)

    inner_surface =
