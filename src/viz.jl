function draw_axes(ax::Axis)
	vlines!(ax, 0.0, color=(:gray40, 0.6), linewidth=1)
	hlines!(ax, 0.0, color=(:gray40, 0.6), linewidth=1)
end

@doc raw"""
    viz_response(data, 
                 title="", xlabel="time, t", 
                 ylabel="output, y(t)",
                 savename=nothing)

plot `data[:, :output]` vs. `data[:, :t]` to visualize the response of a system to an input. 
typically the data frame, `data`, is returned from [`simulate`](@ref).

# Arguments
* `data::DataFrame`: data frame of time series data, containing a `:t` column for times and `:output` column for the outputs.
* `title::String`: title of plot
* `xlabel::String`: x-label
* `ylabel::String`: y-label
* `savename::Union{Nothing, String}`: filename to save as a figure in .png format.

# Returns
`CairoMakie.jl` `Figure` object. this will display in a Pluto.jl notebook.

`CairoMakie.jl` commands can be invoked after `viz_response` to make further changes to the figure panel by e.g.:
```julia
fig = viz_response(data)
ax = current_axis(fig)
ax.xlabel = "new xlabel"
xlims!(ax, 0, 15)
```


# Example
```
g = 4 / (4 * s ^ 2 + 0.8 * s + 1)
u = 1 / s
data = simulate(g * u, (0.0, 50.0))
fig = viz_response(data)
```
"""
function viz_response(data::DataFrame;
                      title::String="",
                      xlabel::String="time, t",
                      ylabel::String="output, y(t)",
                      savename::Union{Nothing, String}=nothing
    )

    fig = Figure()
    ax  = Axis(fig[1, 1], xlabel=xlabel, ylabel=ylabel, title=title)
    draw_axes(ax)
    lines!(data[:, :t], data[:, :output])
    if ! isnothing(savename)
        save(savename, fig, px_per_unit=1)
    end
    return fig
end

"""
    viz_poles_and_zeros(g, savename=nothing)

plot the zeros and poles of the transfer function `g` in the complex plane.

returns a `CairoMakie.jl` `Figure` object for further modification.
"""
function viz_poles_and_zeros(tf::TransferFunction; savename::Union{Nothing, String}=nothing)
    z, p, k = zeros_poles_gain(tf)
    
    fig = Figure()
	ax  = Axis(fig[1, 1], xlabel="Re", ylabel="Im", title="poles and zeros")
    draw_axes(ax)
	scatter!(real.(z), imag.(z), marker=:o, label="zeros", markersize=15)
	scatter!(real.(p), imag.(p), marker=:x, label="poles", markersize=15)
	axislegend()
    if ! isnothing(savename)
        save(savename, fig, px_per_unit=1)
    end
	return fig
end

"""
    nyquist_diagram(tf, nb_pts=500, ω_max=10.0, savename=nothing)

plot the Nyquist diagram for a transfer function `tf` to visualize its frequency response. `s=-1` is plotted as a red `+`. `nb_pts` changes the resolution. `ω_max` gives maximum frequency considered.
"""
function nyquist_diagram(tf::TransferFunction; nb_pts::Int=500, ω_max::Float64=10.0, savename::Union{Nothing, String}=nothing)
    ω_neg = range(-ω_max, 0.0, length=nb_pts)
    ω_pos = range(0.0, ω_max, length=nb_pts)

    g_iω_neg = [evaluate(tf, ω_i * im) for ω_i in ω_neg]
    g_iω_pos = [evaluate(tf, ω_i * im) for ω_i in ω_pos]

    fig = Figure()
	ax  = Axis(fig[1, 1], xlabel="Re[G(iω)]", ylabel="Im[G(iω)]", title="Nyquist diagram")
    draw_axes(ax)
    lines!(real(g_iω_neg), imag(g_iω_neg))
    lines!(real(g_iω_pos), imag(g_iω_pos))
    # plot -1
    scatter([-1], [0], marker=:+, markersize=15, color="red")
    if ! isnothing(savename)
        save(savename, fig, px_per_unit=1)
    end
    return fig
end

"""
    root_locus(g_ol, max_mag_Kc=10.0, nb_pts=500, savename=nothing)

visualize the root locus plot of an open-loop transfer function `g_ol`.

# Arguments
* `g_ol::TransferFunction`: the open-loop transfer function of the closed loop system
* `max_mag_Kc::Float64=10.0`: the maximum magnitude by which the gain of `g_ol` is 
    scaled in order to see the roots traversing the plane
* `nb_pts::Int=500`: the number of gains to explore. increase for higher resolution.
"""
function root_locus(g_ol::TransferFunction;
        max_mag_Kc::Float64=10.0, nb_pts::Int=500, savename=nothing)
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

    fig = Figure()
	ax  = Axis(fig[1, 1], xlabel="Re", ylabel="Im", title="root locus")
    draw_axes(ax)
    # plot poles; corresponds to Kc = 0
    scatter!(real.(p), imag.(p), marker=:x, label="poles", markersize=15, color="black")
    # plot zeros; corresponds to |Kc| → ∞
    if length(z) > 0
        scatter!(real.(z), imag.(z), marker=:o, label="zeros", markersize=15, color="black")
    end
    # plot roots traversing plane
    for i = 1:length(p)
        lines!(real.(rloc[:, i]), imag.(rloc[:, i]))
    end
    axislegend()
    if ! isnothing(savename)
        save(savename, fig, px_per_unit=1)
    end
    return fig
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
    
    fig = Figure()
    axs = [Axis(fig[1, 1], ylabel="|g(iω)|", title="Bode plot"),
           Axis(fig[2, 1], xlabel="ω")]
    linkxaxes!(axs...)
    axs[1].lines!(ω, abs.(g_iω))
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
    tight_layout()
    return axs
end

@doc raw"""
    mk_gif(data, plot_title="", plot_xlabel="time, t", 
                 plot_ylabel="output, y(t)",
                 savename="response")

make a .gif of the process response.
`data` is a data frame with two columns, `:t` and `:output`, likely returned from [`simulate`](@ref).
accepts same arguments as [`viz_response`](@ref).
ImageMagick must be installed to create the .gif.
the .gif is saved as a file `savename`.

# Arguments
* `data::DataFrame`: data frame of time series data, containing a `:t` column for times and `:output` column for the outputs.
* `plot_title::String`: title of plot
* `plot_xlabel::String`: x-label
* `plot_ylabel::String`: y-label
* `savename::String`: filename to save as a .gif. .gif extension automatically appended if not provided.
"""
function mk_gif(t::Array{Float64}, y::Array{Float64};
        plot_title::String="",
        plot_xlabel::String=L"time, $t$",
        plot_ylabel::String=L"output, $y(t)$",
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
