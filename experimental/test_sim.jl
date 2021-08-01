
n0 = 1e13
T0 = 20
beta0 = 0.3
B0 = 22e3
q = 3
Nx = 100
Nz = 1

a = 67
R0 = 167
k = 0.2
dt = 0.001

wrk = rectangle_workspace(Nx, Nz, q, a, R0, k, n0, T0, B0, dt)

lnn = log.(hcat([f_to_grid((x,y)->2*exp(-y^2/(2*0.7^2)), wrk.GRID) for _=1:3]...))
lnTe = log.(hcat([f_to_grid((x,y)->2*exp(-y^2/(2*0.7^2)), wrk.GRID) for _=1:3]...))
lnTi = log.(hcat([f_to_grid((x,y)->2*exp(-y^2/(2*0.7^2)), wrk.GRID) for _=1:3]...))
u = zeros(Float64, (wrk.GRID.Nk*wrk.GRID.Nz, 3))
w = zeros(Float64, (wrk.GRID.Nk*wrk.GRID.Nz, 3))
A = zeros(Float64, (wrk.GRID.Nk*wrk.GRID.Nz, 3))

n = exp.(lnn[:,3])
Te = exp.(lnTe[:,3])
Ti = exp.(lnTi[:,3])
Pe = n .* Te
Pi = n .* Ti
j = zeros(Float64, wrk.GRID.Nk*wrk.GRID.Nz)
jn = j ./ n

Pi_xx = wrk.Dxx * Pe
Pi_yy = wrk.Dyy * Pe
lnn_x = wrk.Dx * lnn[:,2]
lnn_y = wrk.Dy * lnn[:,2]

phi_b = zeros(Float64, wrk.GRID.Nk*wrk.GRID.Nz)
phi = vorticity_eqn(zeros(Float64, wrk.GRID.Nk*wrk.GRID.Nz), w[:,2], Pi_xx, Pi_yy, Te, n, lnn_x, lnn_y, phi_b, wrk.params[2], wrk)
psi = zeros(Float64, wrk.GRID.Nk*wrk.GRID.Nz)

leapfrog!(lnn, lnTe, lnTi, u, w, A, phi, psi, n, Te, Ti, Pe, Pi, j, jn, wrk)
