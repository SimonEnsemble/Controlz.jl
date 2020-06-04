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
g = 4 / (4 * s ^ 2 + 0.8 * s + 1)
u = 1 / s
t, y = simulate(g * u, (0.0, 50.0))
viz_response(t, y)
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
    if ! isnothing(savename)
        tight_layout()
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
    nyquist_diagram(tf, nb_pts=500, ω_max=10.0)

plot the Nyquist diagram for a transfer function `tf` to visualize its frequency response. `s=-1` is plotted as a red `+`. `nb_pts` changes the resolution. `ω_max` gives maximum frequency considered.
"""
function nyquist_diagram(tf::TransferFunction; nb_pts::Int=500, ω_max::Float64=10.0)
    ω_neg = range(-ω_max, 0.0, length=nb_pts)
    ω_pos = range(0.0, ω_max, length=nb_pts)

    g_iω_neg = [evaluate(tf, ω_i * im) for ω_i in ω_neg]
    g_iω_pos = [evaluate(tf, ω_i * im) for ω_i in ω_pos]

    figure()
    plot(real(g_iω_neg), imag(g_iω_neg), zorder=100)
    plot(real(g_iω_pos), imag(g_iω_pos), zorder=100)
    draw_axes()
    # plot -1
    plot([-1], [0], marker="+", color="r", zorder=1000, markersize=15)
    xlabel("Re[G(iω)]")
    ylabel("Im[G(iω)]")
    title("Nyquist diagram")
    tight_layout()
end

"""
    root_locus(g_ol, max_mag_Kc=10.0, nb_pts=500)

visualize the root locus plot of an open-loop transfer function `g_ol`.

# Arguments
* `g_ol::TransferFunction`: the open-loop transfer function of the closed loop system
* `max_mag_Kc::Float64=10.0`: the maximum magnitude by which the gain of `g_ol` is 
    scaled in order to see the roots traversing the plane
* `nb_pts::Int=500`: the number of gains to explore. increase for higher resolution.
"""
function root_locus(g_ol::TransferFunction;
        max_mag_Kc::Float64=10.0, nb_pts::Int=500)
    # compute zeros, poles, and gain of open loop transfer function
    z, p, k = zeros_poles_k(g_ol)

    # decide set of controller gains to explore
    #   make sure Kc is the same sign as gain of G_OL
    #   so G_OL has positive gain
    Kcs = k * 10.0 .^ range(-6, log(max_mag_Kc), length=nb_pts)
    pushfirst!(Kcs, 0.0)

    # store roots of characteristic eqn. here as Kc varies
    rloc = zeros(Complex, length(Kcs), length(p))
    fill!(rloc, NaN)

    # roots of 1 + G_OL when Kc = 0 are poles of G_OL
    rloc[1, :] = p

    # loop thru controller gains, compute roots of 1 + G_OL(s) = 0
    for i = 2:length(Kcs)
        # compute characteristic polynomial 1 + Kc * G_OL(s)
        c_poly = characteristic_polynomial(Kcs[i] * g_ol)
        # compute roots of characteristic polynomial
        roots_c_poly = roots(c_poly)
        # store roots
        entry_filled = [false for _ = 1:length(p)]
        for i_p = 1:length(roots_c_poly)
            ### find closest previous root
            # first, compute distances from each
            distance_to_previous_roots = abs.(rloc[i-1, :] .- roots_c_poly[i_p])
            distance_to_previous_roots[entry_filled] .= Inf
            # get closest
            id_closest_root = argmin(distance_to_previous_roots)
            entry_filled[id_closest_root] = true
            # store root in this entry
            rloc[i, id_closest_root] = roots_c_poly[i_p]
        end
    end

    figure()
    # plot poles; corresponds to Kc = 0
    scatter(real.(p), imag.(p), marker="x", label="poles",
            color="k", s=50, zorder=100)
    # plot zeros; corresponds to |Kc| → ∞
    if length(z) > 0
        scatter(real.(z), imag.(z), marker="o", label="poles",
                color="k", s=50, zorder=100, facecolor="None")
    end
    # plot roots traversing plane
    for i = 1:length(p)
        plot(real.(rloc[:, i]), imag.(rloc[:, i]),
            zorder=100, color="C$(i-1)")
    end
    xlabel("Re")
    ylabel("Im")
    draw_axes()
    title("root locus")
    tight_layout()
end

"""
    axs = bode_plot(tf, log10_ω_min=-4.0, log10_ω_max=4.0, nb_pts=300)

draw the Bode plot of a transfer function `tf` to visualize its frequency response.
returns the two axes of the plot for further tuning via `matplotlib` commands.

adjust the range of frequencies that the Bode plot presents with `log10_ω_min` and `log10_ω_max`.
increase the resolution of the Bode plot with `nb_pts`.
"""
function bode_plot(g::TransferFunction; log10_ω_min::Float64=-3.0, log10_ω_max::Float64=3.0, nb_pts::Int=300)
    ω = 10.0 .^ range(log10_ω_min, log10_ω_max, length=nb_pts)
    g_iω = [evaluate(g, im * ω_i) for ω_i in ω]
    ∠g_iω = zeros(length(g_iω))

    circle_counter = 0
    ∠g_iω[1] = angle(g_iω[1])
    for i = 2:length(g_iω)
        ∠g_iω[i] = angle(g_iω[i]) - circle_counter * 2 * π
        if ∠g_iω[i] - ∠g_iω[i-1] > π
            ∠g_iω[i] -= 2 * π
            circle_counter += 1
        end
    end

    fig, axs = subplots(2, 1, sharex=true, figsize=(8, 7))
    axs[1].plot(ω, abs.(g_iω), color="C1")
    axs[1].set_ylabel(L"$|g(i\omega)|$")
    axs[1].set_title("Bode plot")
    axs[1].set_xscale("log")
    axs[1].set_yscale("log")
    axs[2].set_xscale("log")
    axs[2].plot(ω, ∠g_iω / π, color="C1")
    axs[2].yaxis.set_major_formatter(PyPlot.matplotlib.ticker.FormatStrFormatter(L"%g$\pi$"))
    axs[2].set_ylabel(L"$\angle g(i\omega)$")
    for ax in axs
        ax.minorticks_on()
        ax.grid(b=true, which="minor", alpha=0.25)
    end
    xlabel(L"frequency, $\omega$")
    tight_layout()
    return axs
end

@doc raw"""
    mk_gif(t, y, plot_title="", plot_xlabel="time, t", 
                 plot_ylabel="output, y(t)",
                 savename="response")

make a .gif of the process response.
`t` and `y` are outputs of [`simulate`](@ref).
accepts same arguments as [`viz_response`](@ref).
ImageMagick must be installed to create the .gif.
the .gif is saved as a file `savename`.

# Arguments
* `t::Array{Float64}`: array of times
* `y::Array{Float64}`: array of values of response variables at the corresponding times in `t`
* `plot_title::Union{String, LaTeXString}`: title of plot
* `plot_xlabel::Union{String, LaTeXString}`: x-label
* `plot_ylabel::Union{String, LaTeXString}`: y-label
* `savename::String`: filename to save as a .gif. .gif extension automatically appended if not provided.
"""
function mk_gif(t::Array{Float64}, y::Array{Float64};
        plot_title::Union{String, LaTeXString}="",
        plot_xlabel::Union{String, LaTeXString}=L"time, $t$",
        plot_ylabel::Union{String, LaTeXString}=L"output, $y(t)$",
        savename::String="response.gif"
    )
    if length(t) > 999
        error("too many points and thus images; reduce the number of points")
    end

    # let matplotlib determine x, y lims:
    plot(t, y)
    xmin, xmax, ymin, ymax = axis()
    close()
    
    # name to save image i as.
    step_to_image(i::Int) = @sprintf("__y_of_t_snapshot_%03d.png", i)

    # save series of images
    for i = 2:length(t)
        viz_response(t[1:i], y[1:i], plot_title=plot_title, plot_xlabel=plot_xlabel,
            plot_ylabel=plot_ylabel)
        xlim([xmin, xmax])
        ylim([ymin, ymax])
        tight_layout()
        savefig(step_to_image(i), format="png", dpi=100)
        close()
    end
    
    if ! occursin(".gif", savename)
        savename *= ".gif"
    end
    try
        run(`convert -delay 20 -loop 0 __y_of_t_snapshot_\*.png $savename`)
        @info "see " * savename
    catch
        @warn "You must install ImageMagick for `convert` to run. Exiting..."
    end

    # clean up
    for i = 2:length(t)
        rm(step_to_image(i))
    end
end
