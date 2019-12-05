"""
Convert a polynomial to a string for printing.
"""
function poly_to_string(poly::Poly)
    io = IOBuffer();
    printpoly(io, poly, descending_powers=true)
    return String(take!(io))
end

function Base.show(io::IO, tf::TransferFunction)
    top = poly_to_string(tf.numerator)
    bottom = poly_to_string(tf.denominator)
    
    nb_dashes = maximum([length(top), length(bottom)])
    
    # print numerator
    nb_slack = floor(Int, (nb_dashes - length(top)) / 2) # for centering
    for i = 1:nb_slack
        print(io, " ")
    end
    println(io, top)
    
    # print dashes to separate numerator and denominator
    for _ = 1:nb_dashes
        print(io, "-")
    end
    if tf.time_delay != 0.0
        print(io, " e^(-$(tf.time_delay)*s)")
    end
    print(io, "\n")
    
    # print denominator
    nb_slack = floor(Int, (nb_dashes - length(bottom)) / 2) # for centering
    for i = 1:nb_slack
        print(io, " ")
    end
    println(io, bottom)
end
