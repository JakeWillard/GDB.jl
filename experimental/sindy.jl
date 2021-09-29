
function sparse_regression(X, Theta, lambda, n)

    N = size(X)[2]

    # initial guess
    C = Theta \ X

    # find coefficients and apply threshold
    i = findall(abs.(C) .< lambda)
    j = findall(abs.(C) > lambda)
    C[i] .= 0

    for _=1:n

        for k=1:N
            C[j[:,k],k] = Theta[:,j[:,k]] \ X[:,k]
        end

        # find coefficients and apply threshold
        i = findall(abs.(C) .< lambda)
        j = findall(abs.(C) > lambda)
        C[i] .= 0
    end

    return C
end
