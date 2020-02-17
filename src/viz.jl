using LaTeXStrings

function draw_axes()
    axvline(x=0, color="0.8", lw=2, zorder=1)
    axhline(y=0, color="0.8", lw=2, zorder=1)
end

"""
    viz_response(t, y, 
                 plot_title="", plot_xlabel="time, t", 
                 plot_ylabel="output, y(t)",
                 savename=nothing)

plot `y` vs. `t` to visualize the response of a system to an input. typically `t` and `y` are outputs of [`simulate`](@ref).

note that PyPlot.jl (matplotlib) commands can be invoked after `viz_response` to make further changes to the figure panel.
e.g. `xlim([0, 1])` can be applied after `viz_response`.

# Arguments
* `t::Array{Float64}`: array of times
* `y::Array{Float64}`: array of values of response variables at the corresponding times in `t`
* `plot_title::Union{String, LaTeXString}`: title of plot
* `plot_xlabel::Union{String, LaTeXString}`: x-label
* `plot_ylabel::Union{String, LaTeXString}`: y-label
* `savename::Union{Nothing, String}`: filename to save as a figure in .png format (dpi 250).

# Example
```
julia> g = 4 / (4 * s ^ 2 + 0.8 * s + 1)
julia> u = 1 / s
julia> t, y = simulate(g * u, (0.0, 50.0))
julia> viz_response(t, y)
```
"""
function viz_response(t::Array{Float64}, y::Array{Float64}; 
        plot_title::Union{String, LaTeXString}="", 
        plot_xlabel::Union{String, LaTeXString}=L"time, $t$",
        plot_ylabel::Union{String, LaTeXString}=L"output, $y(t)$",
        savename::Union{Nothing, String}=nothing
    )
    
    figure()
    plot(t, y, zorder=100)
    xlabel(plot_xlabel)
    ylabel(plot_ylabel)
    title(plot_title)
    draw_axes()
    tight_layout()
    if ! isnothing(savename)
        if ! occursin(".png", savename)
            savename *= ".png"
        end
        savefig(savename, dpi=250, format="png")
    end
    return nothing
end

"""
    viz_poles_and_zeros(tf)

plot the zeros and poles of the transfer function `tf` in the complex plane.
"""
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
    tight_layout()
    return nothing
end

"""
    nyquist_diagram(tf)

plot the Nyquist diagram for a transfer function `tf` to visualize its frequency response.
"""
function nyquist_diagram(tf::TransferFunction; nb_pts::Int=300)
    ω = range(-10.0, 10.0, length=nb_pts)

    g_iω = [evaluate(tf, ω_i * im) for ω_i in ω]

    figure()
    plot(real(g_iω), imag(g_iω), zorder=100)
    draw_axes()
    xlabel("Re[G(iω)]")
    ylabel("Im[G(iω)]")
    title("Nyquist diagram")
    tight_layout()
    return nothing
end

"""
    root_locus(g_ol)

visualize the root locus plot of an open-loop transfer function `g_ol`.
"""
function root_locus(g_ol::TransferFunction)
    z, p, k = zeros_poles_k(g_ol)

    # keep gain of G_ol(s) +ve, scale with gain K
    Kcs = k * 10.0 .^ range(-6, 2, length=500)
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
    xlabel("Re")
    ylabel("Im")
    draw_axes()
    title("root locus")
    tight_layout()
    return nothing
end

"""
    bode_plot(tf, log10_ω_min=-4.0, log10_ω_max=4.0)

draw the Bode plot of a transfer function `tf` to visualize its frequency response.
"""
function bode_plot(g::TransferFunction; log10_ω_min::Float64=-4.0, log10_ω_max::Float64=4.0)
    ω = 10.0 .^ range(log10_ω_min, log10_ω_max, length=300)
    g_iω = [evaluate(g, im * ω_i) for ω_i in ω]

    fig, (ax1, ax2) = subplots(2, 1, sharex=true, figsize=(8, 7))
    ax1.plot(ω, abs.(g_iω))
    ax1.set_ylabel(L"$|g(i\omega)|$")
    ax1.set_title("Bode plot")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax2.set_xscale("log")
    ax2.plot(ω, angle.(g_iω), color="C1")
    ax2.set_ylabel(L"$\angle g(i\omega)$")
    xlabel(L"frequency, $\omega$")
    tight_layout()
    return nothing
end

"""
    mk_gif(t, y)

make a .gif of the process response. 
`t` and `y` are outputs of [`simulate`](@ref).
accepts same arguments as [`viz_response`](@ref).
ImageMagick must be installed to create the .gif.
"""
function mk_gif(t::Array{Float64}, y::Array{Float64};
        plot_title::Union{String, LaTeXString}="",
        plot_xlabel::Union{String, LaTeXString}=L"time, $t$",
        plot_ylabel::Union{String, LaTeXString}=L"output, $y(t)$",
        savename::Union{Nothing, String}=nothing
    )
    if length(t) > 999
        error("too many points and thus images; reduce the number of points")
    end

    # let matplotlib determine x, y lims:
    plot(t, y)
    xmin, xmax, ymin, ymax = axis()
    close()

    step_to_image(i::Int) = @sprintf("__y_of_t_snapshot_%03d.png", i)

    # same series of images
    for i = 2:length(t)
        viz_response(t[1:i], y[1:i], plot_title=plot_title, plot_xlabel=plot_xlabel,
            plot_ylabel=plot_ylabel)
        xlim([xmin, xmax])
        ylim([ymin, ymax])
        savefig(step_to_image(i), format="png", dpi=100)
        close()
    end

    if isnothing(savename)
        savename = "reponse"
    end
    try
        run(`convert -delay 20 -loop 0 __y_of_t_snapshot_*.png $savename.gif`)
        @info "see " * savename * ".gif"
    catch
        @warning "You must install ImageMagick for `convert` to run. Exiting..."
    end

    # clean up
    for i = 2:length(t)
        rm(step_to_image(i))
    end
end
