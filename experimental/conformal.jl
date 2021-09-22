
function circle_to_shape(z, alpha, kappa)

    _a = 0.5*(z + 1/z)
    _b = 0.5*(z - 1/z)
    return _a*cosh(alpha*_b) + _b*sinh(alpha*_b) + kappa*_b
end


function diff_circle_to_shape(z, alpha, kappa)

    _a = cosh(alpha*(z - 1/z)/2)
    _b = sinh(alpha*(z - 1/z)/2)
    _c = _a + _b
    _d = _a - _b
    _e = (_c + kappa) - (_d - kappa)/z^2
    _f = 0.5*alpha*(z + 1/z)*_c
    _g = 0.5*alpha*(1/z + 1/z^3)*_d

    return (_f + _g) / 2
end


function shape_to_circle(w, alpha, kappa)

    r = circle_to_shape(w, alpha, kappa) - w
    z = w

    xs = zeros(10)
    ys = zeros(10)
    for i=1:10
        z -= r / diff_circle_to_shape(z, alpha, kappa)
        r = circle_to_shape(z, alpha, kappa) - w
        xs[i] = real(z)
        ys[i] = imag(z)
    end
    # while abs(r) > 1e-8
    #     z -= r / diff_circle_to_shape(z, alpha, kappa)
    #     r = circle_to_shape(z, alpha, kappa) - w
    # end

    p = plot(xs, ys)
    return p
end
