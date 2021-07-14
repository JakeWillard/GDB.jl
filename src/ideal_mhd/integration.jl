
@variable "LogDensity"
@variable "Flux"
@variable "Stream"
@variable "Vorticity"
@variable "Current"

LogDensity(x::Vector{Float64}) = LogDensity{Float64, 2}(hcat(x, x))
Flux(x::Vector{Float64}) = Flux{Float64, 2}(hcat(x, x))
Stream(x::Vector{Float64}) = Stream{Float64, 2}(hcat(x, x))
Vorticity(x::Vector{Float64}) = Vorticity{Float64, 2}(hcat(x, x))
Current(x::Vector{Float64}) = Current{Float64, 2}(hcat(x, x))


"""

GOVERNING EQUATIONS (REDUCED MHD IN 2 DIMENSIONS):

lnn_t = -{phi, lnn_t}
w_t = exp(-lnn) * (-{phi, w} + {psi, j})
psi_t = -{phi, psi} + eta*j
j = del^2 * psi
w = del^2 * phi

"""


function partial_t(lnn::LogDensity, phi::Stream, stp::Setup)

    # compute poisson bracket
    lnn_t = (stp.Dx*lnn[:,2]) .* (stp.Dy*phi[:,2]) .- (stp.Dx*phi[:,2]) .* (stp.Dy*lnn[:,2])

    # density source on x-flux surface
    lnn_src = 0.3
    lnn_t = stp.P1 * lnn_t[:] + lnn_src*diag(I - stp.P1)

    # density sink on y-flux surface
    lnn_t = stp.P2 * lnn_t[:]

    # Neumann conditions on hard boundaries
    lnn_t = stp.P3 * lnn_t[:] + 0.5*(I - stp.P3)*(stp.R3 - I)*lnn[:,2]

    return lnn_t
end


function partial_t(w::Vorticity, phi::Stream, psi::Flux, j::Current, stp::Setup)

    # compute poisson bracket
    w_t = (stp.Dx*w[:,2]) .* (stp.Dy*phi[:,2]) .- (stp.Dx*phi[:,2]) .* (stp.Dy*w[:,2])

    # compute current term
    w_t += (stp.Dx*psi[:,2]) .* (stp.Dy*j[:,2]) .- (stp.Dx*j[:,2]) .* (stp.Dy*psi[:,2])

    # Dirichlet conditions on every surface
    w_t = stp.P1*w_t + 0.5*(I - P1)*(R1 + I)*w[:,2]
    w_t = stp.P2*w_t + 0.5*(I - P1)*(R2 + I)*w[:,2]
    w_t = stp.P3*w_t + 0.5*(I - P1)*(R3 + I)*w[:,2]

    return w_t
end


function partial_t(psi::Flux, phi::Stream, j::Current, eta::Float64, stp::Setup)

    # compute poisson bracket
    psi_t = (stp.Dx*psi[:,2]) .* (stp.Dy*phi[:,2]) .- (stp.Dx*phi[:,2]) .* (stp.Dy*psi[:,2])

    # compute current term
    psi_t += eta * j[:,2]

    # Dirichlet conditions on every surface
    psi_t = stp.P1*psi_t + 0.5*(I - P1)*(R1 + I)*psi[:,2]
    psi_t = stp.P2*psi_t + 0.5*(I - P1)*(R2 + I)*psi[:,2]
    psi_t = stp.P3*psi_t + 0.5*(I - P1)*(R3 + I)*psi[:,2]

    return psi_t
end


function solve_vorticity_eqn(phi_b::Stream, w::Vorticity, stp::Setup)

    # clamp to phi_b (given boundary value) on flux surfaces
    A = stp.P1 * stp.L + 0.5*(I - stp.P1)*(stp.R1 + I)
    A = stp.P2 * A + 0.5*(I - stp.P2)*(stp.R2 + I)
    rhs = stp.P1*w[:,2] + 0.5*(I - stp.P1)*(stp.R1 + I)*phi_b[:]
    rhs = stp.P2*rhs[:] + 0.5*(I - stp.P2)*(stp.R2 + I)*phi_b[:]

    # no-slip conditions on hard surfaces
    A = stp.P3 * A + 0.5*(I - stp.P3)*(stp.R3 - I)
    rhs = stp.P3*rhs[:]

    # turn into least-squares problem
    rhs = transpose(A) * rhs[:]
    A = transpose(A) * A

    # solve with backslash
    return A \ rhs
end


function solve_diffusion(f::Union{Vorticity, Flux}, Cdiff::Float64, stp::Setup)

    A = stp.P1*(I - Cdiff*stp.L) + 0.5*(I - P1)*(R1 + I)
    A = stp.P2*A + 0.5*(I - P2)*(R2 + I)
    A = stp.P3*A + 0.5*(I - P3)*(R3 + I)
    rhs = (stp.P1*stp.P2*stp.P3) * f[:,2]

    # turn into least-squares problem
    rhs = transpose(A) * rhs[:]
    A = transpose(A) * A

    # solve with backslash
    return A \ rhs
end


function solve_diffusion(lnn::LogDensity, Cdiff::Float64, stp::Setup)

    A = stp.P1*(I - Cdiff*stp.L) + 0.5*(I - P1)*(R1 + I)
    A = stp.P2*A + 0.5*(I - P2)*(R2 + I)
    A = stp.P3*A + 0.5*(I - P3)*(R3 + I)

    lnn_src = 0.3
    rhs = stp.P1 * lnn[:,2] + lnn_src*diag(I - stp.P1)
    rhs = (stp.P2*stp.P3) * rhs[:]
end


function leapfrog!(lnn::LogDensity, w::Vorticity, psi::Flux, phi::Stream, j::Current, phi_b::Stream, dt::Float64, eta::Float64, Cdiff::Float64, stp::Setup)

    # compute time derivatives
    lnn_t = partial_t(lnn, phi, stp)
    w_t = partial_t(w, phi, psi, j, stp)
    psi_t = partial_t(psi, phi, j, eta, stp)

    # compute predictor step
    lnn_p = LogDensity(lnn[1,:] + dt*lnn_t)
    w_p = Vorticity(w[1,:] + dt*w_t)
    psi_p = Flux(psi[1,:] + dt*psi_t)
    phi_p = Stream(solve_vorticity_eqn(phi_b, w_p, stp))
    j_p = Current(stp.L * psi_p)

    # recompute time derivatives
    lnn_t_new = partial_t(lnn_p, phi_p, stp)
    w_t_new = partial_t(w_p, phi_p, psi_p, j_p, stp)
    psi_t_new = partial_t(psi_p, phi_p, j_p, eta, stp)

    # compute final values
    lnn_f = lnn[2,:] + dt*lnn_t_new
    w_f = w[2,:] + dt*w_t_new
    psi_f = psi[2,:] + dt*psi_t_new
    lnn[1,:] = lnn[2,:]
    lnn[2,:] = lnn_f[:]
    w[1,:] = w[2,:]
    w[2,:] = w_f[:]
    psi[1,:] = psi[2,:]
    psi[2,:] = psi_f[:]
    phi[1,:] = phi[2,:]
    phi[2,:] = solve_vorticity_eqn(phi_b, w, stp)
    j[1,:] = j[2,:]
    j[2,:] = stp.L * psi[2,:]
end
