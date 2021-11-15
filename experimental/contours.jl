

module ContourWrapper

    export contoursw, levels, lines, coordinates

    using Contour

    function contoursw(args...)
        contours(args...)
    end
end
