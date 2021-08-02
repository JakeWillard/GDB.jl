
struct Assets

    GRID :: Grid
    Dx :: SparseMatrixCSC
    Dy :: SparseMatrixCSC
    Dxy :: SparseMatrixCSC
    Dxx :: SparseMatrixCSC
    Dyy :: SparseMatrixCSC
    Dxxx :: SparseMatrixCSC
    Dyyy :: SparseMatrixCSC
    Dxxy :: SparseMatrixCSC
    Dxyy :: SparseMatrixCSC
    Ds :: SparseMatrixCSC
    Dss :: SparseMatrixCSC
    DIFF_lnn :: LinearLeftHandSide
    DIFF_lnTe :: LinearLeftHandSide
    DIFF_lnTi :: LinearLeftHandSide
    DIFF_u :: LinearLeftHandSide
    DIFF_w :: LinearLeftHandSide
    DIFF_A :: LinearLeftHandSide
    HHOLTZ :: LinearLeftHandSide
    P0 :: SparseMatrixCSC
    P1 :: SparseMatrixCSC
    P2 :: SparseMatrixCSC
    P3 :: SparseMatrixCSC
    R1 :: SparseMatrixCSC
    R2 :: SparseMatrixCSC
    R3 :: SparseMatrixCSC
    LAM :: SparseMatrixCSC
    DCHLT1 :: SparseMatrixCSC
    DCHLT2 :: SparseMatrixCSC
    DCHLT3 :: SparseMatrixCSC
    NMANN1 :: SparseMatrixCSC
    NMANN2 :: SparseMatrixCSC
    NMANN3 :: SparseMatrixCSC
    FLXAVG :: SparseMatrixCSC
    TRGT :: Vector{Float64}
    Sn :: Vector{Float64}
    STe :: Vector{Float64}
    STi :: Vector{Float64}
    params :: Vector{Float64}
    dt :: Float64
    N_subcycle :: Int64
    seed :: Vector{Float64}

end
