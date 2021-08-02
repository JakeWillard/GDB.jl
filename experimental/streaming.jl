struct StreamTensor

    E :: SparseMatrixCSC
    SE :: SparseMatrixCSC
    S :: SparseMatrixCSC
    SW :: SparseMatrixCSC
    W :: SparseMatrixCSC
    NW :: SparseMatrixCSC
    N :: SparseMatrixCSC
    NE :: SparseMatrixCSC

end


function StreamTensor(Nx, Ny)

    Fx = spdiagm(-1 => ones(Nx - 1))
    Fx[1,Nx] = 1
    Bx = spdiagm(1 => ones(Nx - 1))
    Bx[Nx,1] = 1

    Fy = spdiagm(-1 => ones(Ny - 1))
    Fy[1,Ny] = 1
    By = spdiagm(1 => ones(Ny - 1))
    By[Ny,1] = 1

    Ix = sparse(I, Nx, Nx)
    Iy = sparse(I, Ny, Ny)

    E = kron(Iy, Fx)
    SE = kron(By, Fx)
    S = kron(By, Ix)
    SW = kron(By, Bx)
    W = kron(Iy, Bx)
    NW = kron(Fy, Bx)
    N = kron(Fy, Ix)
    NE = kron(Fy, Fx)

    return StreamTensor(E, SE, S, SW, W, NW, N, NE)
end


function Base.:*(S::StreamTensor, f::Matrix)

    out = zeros(size(f))

    out[:,1] = S.N * f[:,1]
    out[:,2] = S.NE * f[:,2]
    out[:,3] = S.E * f[:,3]
    .
    .
    .
    .

    return out
end
