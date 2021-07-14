
function relax(M, x, f, i, j, n)

    out = zeros(Float64, size(x))
    y = zeros(Float64, size(x))
    y[:] = x[:]
    for dummy=1:n-1
        y[i:j] = M[i:j,:] * y + f[i:j]
    end
    # if these are being added eventually, we only want the the portion of the matrix corresponding to the relevant region to be copied
    # with the loop before, this likely performs the operation n times before returning the output
    # out[i:j] = M[i:j,:] * y + f[i:j]

    # my intuition would be
    out[i:j] = M[i:j,:] * x + f[i:j]
    # assuming it is already being computed N times in the below function
    return sparse(out)
end


function p_jacobi(A::SparseMatrixCSC, x::Vector{Float64}, b::Vector{Float64}, w::Float64, N::Int, n::Int, chunks::Int, threshold::Float64)

    # inverse of matrix constructed from diagonal of A
    Dinv = inv(Diagonal(A))
    
    # w = weight factor, w=1 for standard Jacobi iteration
    # defining M as the identity matrix minus the inverse diagonal matrix multiplied by the original matrix
    M = I - w*Dinv*A
    # defining f as the inverse of the matrix D, weighted with w 
    f = w*Dinv * b

    bnorm = norm(b)

    # takes next higher integer value of different values in the range from 1 to the length of x, where chunks gives how many rows you want to assign to a node
    nodes = Int32[ceil(i) for i in LinRange(1, length(x), chunks)]
    # xnext - values of x in sequential iteration
    xnext = zeros(Float64, size(x))



    # predetermined error threshold - iteration ceases while error is below this value
    err = 1.0
    xnext[:] = x[:]
    while err > threshold 
        # iterating over a dummy variable N times - only necessary if you want to iterate N times before checking the error against threshold
        for dummy=1:N
            # serial
            xnext[:] = M * xnext + f # actual equation describing the Jacobi iteration
            # parallel
            xnext[:] = sum(pmap(i -> relax(M, xnext, f, nodes[i], nodes[i+1]-1, n) : i, 1:chunks-1)) # shouldn't use both versions at once
            # how to get rid of repeated indices without getting rid of last element?
        end
        

        # potential way to do asynchronous calculation - calculate each chunk separately
        # way to index by chunk length?
        # done implicitly in pmap function
        """@distributed for i in nodes-1
            xnext[(i-1)*chunks:(i*chunks-1)] = M[(i-1)*chunks:(i*chunks-1),:] * x + f[(i-1)*chunks:(i*chunks-1)]
        end"""
        

        
        err = norm(A * xnext - b) / bnorm # definition of the error function
        # loop iterates N times then checks if error is greater/less than threshold, then does it again if it is still greater than threshold
        x[:] = xnext[:]
    end
    return x


end


function jacobi_smooth(A::SparseMatrixCSC, x::Vector{Float64}, b::Vector{Float64}, w::Float64, N::Int, n::Int, chunks::Int)
    # inverse of matrix constructed from diagonal of A
    Dinv = inv(Diagonal(A))
    
    # w = weight factor, w=1 for standard Jacobi iteration
    # defining M as the identity matrix minus the inverse diagonal matrix multiplied by the original matrix
    M = I - w*Dinv*A
    # defining f as the inverse of the matrix D, weighted with w 
    f = w*Dinv * b

    bnorm = norm(b)

    # takes next higher integer value of different values in the range from 1 to the length of x, where chunks gives how many rows you want to assign to a node
    nodes = Int32[ceil(i) for i in LinRange(1, length(x), chunks)]
    # xnext - values of x in sequential iteration
    xnext = zeros(Float64, size(x))


    xnext[:] = x[:]
    # iterating over a dummy variable N times - only necessary if you want to iterate N times before checking the error against threshold
    for dummy=1:N
        # serial
        xnext[:] = M * xnext + f # actual equation describing the Jacobi iteration
        # parallel
        xnext[:] = sum(pmap(i -> relax(M, xnext, f, nodes[i], nodes[i+1]-1, n) : i, 1:chunks-1)) # shouldn't use both versions at once
        # how to get rid of repeated indices without getting rid of last element?
    end
       
    err = norm(A * xnext - b) / bnorm # definition of the error function
    # loop iterates N times then checks if error is greater/less than threshold, then does it again if it is still greater than threshold
    x[:] = xnext[:]
    return x
end