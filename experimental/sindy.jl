
function sparse_regression(X, Theta, lambda; maxiters=20)

    N = size(X)[2]

    # initial guess
    C = Theta \ X

    # find coefficients and apply threshold
    i = abs.(C) .< lambda
    j = 1 .- i
    C[i] .= 0
    sparsity = sum(i[:])

    for _=1:maxiters

        for k=1:N
            C[j[:,k],k] = Theta[:,j[:,k]] \ X[:,k]
        end

        # find coefficients and apply threshold
        i = abs.(C) .< lambda
        j = 1 .- i
        C[i] .= 0

        # if the sparsity hasn't changed, stop iteration
        new_sparsity = sum(i[:])
        if new_sparsity == sparsity
            break
        end

        sparisty = new_sparsity
    end

    return C
end
