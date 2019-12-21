using LaTeXStrings
using PyPlot
# hipster plot theme
PyPlot.matplotlib.style.use("hipster.mplstyle")

function draw_axes()
    axvline(x=0, color="0.7", lw=2, zorder=1)
    axhline(y=0, color="0.7", lw=2, zorder=1)
end

function viz_response(t::Array{Float64}, y::Array{Float64}; 
        plot_title::Union{String, LaTeXString}="", 
        plot_xlabel::Union{String, LaTeXString}=L"time, $t$",
        plot_ylabel::Union{String, LaTeXString}=L"output, $y(t)$",
    )
    
    figure()
    plot(t, y, zorder=100)
    xlabel(plot_xlabel)
    ylabel(plot_ylabel)
    title(plot_title)
    draw_axes()
    return nothing
end

function viz_poles_and_zeros(tf::TransferFunction)
    z, p, k = zeros_poles_gain(tf)
    
    figure()
    scatter(real.(z), imag.(z), marker="o", label="zeros", color="C1", zorder=100, s=50)
    scatter(real.(p), imag.(p), marker="x", label="poles", color="C2", s=50, zorder=100)
    legend()
    draw_axes()
    xlabel("Re")
    ylabel("Im")
    title("poles and zeros")
    return nothing
end

function viz_nyquist_diagram(tf::TransferFunction; nb_pts::Int=200)
    ω = range(-10.0, 10.0, length=nb_pts)

    g_iω = [evaluate(tf, ω_i * im) for ω_i in ω]

    figure()
    plot(real(g_iω), imag(g_iω), zorder=100)
    draw_axes()
    xlabel("Re[G(iω)]")
    ylabel("Im[G(iω)]")
    title("Nyquist diagram")
    return nothing
end

function viz_root_locus(g_ol::TransferFunction)
    z, p, k = zeros_poles_gain(g_ol)

    # keep gain of G_ol(s) +ve, scale with gain K
    Kcs = k * 10.0 .^ range(-6, 2, length=300)
    pushfirst!(Kcs, 0.0)

    rloc = zeros(Complex, length(Kcs), length(p))

    for (i_k, Kc) in enumerate(Kcs)
        c_poly = characteristic_polynomial(Kc * g_ol)
        roots_c_poly = roots(c_poly)
        for i_p = 1:length(roots_c_poly)
            rloc[i_k, i_p] = roots_c_poly[i_p]
        end
    end

    figure()
    for i = 1:length(p)
        plot(real.(rloc[:, i]), imag.(rloc[:, i]), zorder=100, color="C$(i-1)")
        scatter(real.([p[i]]), imag.([p[i]]), marker="x", label="poles", color="C$(i-1)", s=50, zorder=100)
    end
    xlabel("Re"),
    ylabel("Im"),
    draw_axes()
    title("root locus")
    return nothing
end
