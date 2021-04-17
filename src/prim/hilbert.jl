
# NOTE: This does what its supposed to do but I don't think its
# a good idea to use it actually.
# NOTE: I just copied this from wikipedia, not sure how it actually works
function hilbert_ordering(d::Int, n::Int)

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
