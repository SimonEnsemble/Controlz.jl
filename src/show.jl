"""
Convert a polynomial to a string for printing.
"""
function poly_to_string(poly::Polynomial)
    io = IOBuffer();
    printpoly(io, poly, descending_powers=true)
    return String(take!(io))
end

function show_time_delay(io::IO, time_delay::Union{Int, Float64})
    if time_delay != 0.0
        print(io, " e^($(-time_delay)*s)")
    end
end

function Base.show(io::IO, tf::TransferFunction)
    top = poly_to_string(tf.numerator)
    bottom = poly_to_string(tf.denominator)

    print(io, "\n")

    # is this a rational polynomial? then just print the numerator / constant
    if degree(tf.denominator) == 0
        # normalized numerator
        p = tf.numerator / tf.denominator[0]
        print(io, poly_to_string(p))
        show_time_delay(io, tf.time_delay)
        return nothing
    end
    
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
    show_time_delay(io, tf.time_delay)
    print(io, "\n")

    # print denominator
    nb_slack = floor(Int, (nb_dashes - length(bottom)) / 2) # for centering
    for i = 1:nb_slack
        print(io, " ")
    end
    print(io, bottom)
    return nothing
end

function Base.show(io::IO, pc::PController)
    println(io, "P controller")
    println(io, "\tcontroller gain Kc = ", pc.Kc)
end

function Base.show(io::IO, pic::PIController)
    println(io, "PI controller")
    println(io, "\tcontroller gain Kc = ", pic.Kc)
    println(io, "\tintegral time constant τI = ", pic.τI)
end

function Base.show(io::IO, pidc::PIDController)
    println(io, "PID controller")
    println(io, "\tcontroller gain Kc = ", pidc.Kc)
    println(io, "\tintegral time constant τI = ", pidc.τI)
    println(io, "\tderivative time constant τD = ", pidc.τD)
    println(io, "\tderivative filter α = ", pidc.α)
end

function Base.show(io::IO, m::Margins)
    println(io, "-- gain/phase margin info--")
    @printf(io, "\tcritical frequency ω_c [rad/time]:       %.5f\n", m.ω_c)
    @printf(io, "\tgain crossover frequency ω_g [rad/time]: %.5f\n", m.ω_g)
    @printf(io, "\tgain margin:                             %.5f\n", m.gain_margin)
    @printf(io, "\tphase margin:                            %.5f\n", m.phase_margin)
end
