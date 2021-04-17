

# NOTE: This does what its supposed to do but I don't think its
# a good idea to use it actually.
# NOTE: I just copied this from wikipedia, not sure how it actually works
function hilbert_curve(d::Int, n::Int)

    t = d
    x = 0.0
    y = 0.0

    for i=1:n
        s = 2^(i-1)
        rx = 1 & div(t, 2)
        ry = 1 & xor(t, rx)

        if (ry == 0)
            if (rx == 1)
                x = s - 1 - x
                y = s - 1 - y
            end

            r = Int[y, x]
            x = r[1]
            y = r[2]
        end

        x += s * rx
        y += s * ry
        t = div(t, 4)
    end

    return Int[x, y]
end


struct Grid

    points::Matrix{Float64}
    projection::SparseMatrixCSC
    inverse_projection::SparseMatrixCSC

end


function Grid(inside::Function, n::Int)

    N = 2^n
    Nx = Int(sqrt(N))
    Nk = 0
    points = zeros(Float64, (2, N))

    proj_rows = Int32[i for i=1:N]
    proj_cols = zeros(Int32, N)
    proj_vals = ones(Int32, N)

    for i=1:Nx
        for j=1:Nx
            x = i / Nx
            y = j / Nx
            if inside([x, y])
                Nk += 1
                points[:,Nk] = [x, y]
                proj_cols[Nk] = i + Nx*(j-1)
            end
        end
    end

    P = sparse(proj_rows[1:Nk], proj_cols[1:Nk], proj_vals[1:Nk], Nk, N)
    Pinv = transpose(P)

    return Grid(points, P, Pinv)

end
