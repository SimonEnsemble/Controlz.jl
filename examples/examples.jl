### A Pluto.jl notebook ###
# v0.17.4

using Markdown
using InteractiveUtils

# ╔═╡ 4952d2aa-6b42-11ec-15ea-ed5764c1810b
import Pkg; Pkg.activate()

# ╔═╡ 5da3e954-e940-45a7-b355-45e78b5ce38b
using Controlz, DataFrames, CairoMakie, ColorSchemes, Colors

# ╔═╡ ecab99c5-d929-4a0a-9f81-89ae382aaf45
set_theme!(cool_theme)

# ╔═╡ 2136ca01-87b4-47e9-93b1-25892896413b
md"
## First order plus time delay, response to step

$g(s)=\dfrac{Ke^{-\theta s}}{\tau s + 1}$

response to step input $U(s)=1/s$.
"

# ╔═╡ e664515e-d86b-4a91-9d47-2d77ad19d06a
function sim_first_order_response()
	K = 2.0 # gain
	τ = 4.0 # time constant
	θ = 1.5 # time delay
	
	g = K * exp(-θ * s) / (τ * s + 1) # FOPTD transfer function

	U = 1 / s # step input
	Y = g * U # response

	t_final = 15.0
	data = simulate(Y, t_final)

	fig = viz_response(data, title="FOPTD step response", 
	                   savename="../docs/src/FOPTD_step_response.png")
	return fig
end

# ╔═╡ e2f12441-c8e7-486c-99b4-0ba302d41895
sim_first_order_response()

# ╔═╡ 5bcd4573-7656-4ebe-ba88-45adf2cc4127
md"
## Second order underdamped, response to step

$g(s)=\dfrac{4}{4 s ^2 + 0.8 s + 1}$

response to step input $U(s)=1/s$.
"

# ╔═╡ a25db436-2022-4c33-87fc-36f9b507823f
function sim_second_order_response()
	g = 4 / (4 * s ^ 2 + 0.8 * s + 1)
	
	U = 1 / s
	Y = g * U
	
	data = simulate(Y, 50.0)
	
	fig = viz_response(data, title="SO underdamped step response", 
					   savename="../docs/src/SO_underdamped_step_response.png")
	return fig
end

# ╔═╡ b4433a1c-7d03-440f-af08-51d6d5f6d822
sim_second_order_response()

# ╔═╡ f52a6b52-8202-4c3f-a7b4-4280ffc191d4
md"## viz poles/zeros of a transfer function"

# ╔═╡ 4da06cb3-5086-4fe8-b728-a3f8e7e889c7
function inspect_poles_and_zeros()
	g = (s + 2) / (s^2 + 1/4)
	viz_poles_and_zeros(g, savename="../docs/src/example_poles_and_zeros.png")
end

# ╔═╡ 4cb723a0-d3de-4a98-8d41-0e045a161d07
inspect_poles_and_zeros()

# ╔═╡ 252f5036-62f7-45ef-b7ca-5865aab4cda8
md"## Nyquist diagram"

# ╔═╡ a5c6b75a-be4b-4c82-ab87-2b72463b6303
function draw_nyquist_diagram()
	g = 1 / (s^2 + s + 1) # https://en.wikipedia.org/wiki/Nyquist_stability_criterion
	fig = nyquist_diagram(g, savename="../docs/src/example_nyquist.png")
	return fig
end

# ╔═╡ 068905e6-9cde-4cc9-a042-a36b79c1a234
draw_nyquist_diagram()

# ╔═╡ 8a3d4776-5968-450d-982f-00249caa82ae
md"## root locus diagram"

# ╔═╡ bbd279b4-184d-4e6f-bacf-b468d5c1fd36
function draw_root_locus_diagram()
	g_ol = 4 / (s + 3) / (s + 2) / (s + 1)
	
	fig = root_locus(g_ol, savename="../docs/src/example_root_locus.png")
	return fig
end

# ╔═╡ 24a1c30c-a21c-471d-bcf7-d90fe7200a8c
draw_root_locus_diagram()

# ╔═╡ afee73f3-740d-4114-96c5-b50b1e025e6c
md"## Bode plot"

# ╔═╡ d3996cd9-6890-4445-b2fe-c38ecc76ce77
function draw_bode_plot()
	g = 3 / (s+1)

	fig = bode_plot(g, savename="../docs/src/example_bode.png")
	return fig
end

# ╔═╡ 26d9f8ea-5d85-4cdd-b500-0e2d63f60e81
draw_bode_plot()

# ╔═╡ 13ce095f-2ba2-42a5-b753-cd6a43f024af
function mybode_plot(g::TransferFunction; log10_ω_min::Float64=-3.0, log10_ω_max::Float64=3.0, nb_pts::Int=300)
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

    fig = Figure(resolution=(800, 600))
    axs = [Axis(fig[1, 1], xscale=log10, yscale=log10,
		        ylabel="|g(iω)|", title="Bode plot",
				xminorticksvisible=true, xminorgridvisible=true,
				yminorticksvisible=true, yminorgridvisible=true),
           Axis(fig[2, 1], xscale=log10,
			    xlabel="ω", ylabel="∠g(iω)",
		        xminorticksvisible=true, xminorgridvisible=true,
				yminorticksvisible=true, yminorgridvisible=true)
			]
    # linkxaxes!(axs...)
    lines!(axs[1], ω, abs.(g_iω))
    lines!(axs[2], ω, ∠g_iω / π)
	# axs[2].yaxis.set_major_formatter(PyPlot.matplotlib.ticker.FormatStrFormatter(L"%g$\pi$"))
    return fig
end


# ╔═╡ 7310762a-3c02-42e0-bbb7-85011ba4e144
mybode_plot(3/(s+1))

# ╔═╡ Cell order:
# ╠═4952d2aa-6b42-11ec-15ea-ed5764c1810b
# ╠═5da3e954-e940-45a7-b355-45e78b5ce38b
# ╠═ecab99c5-d929-4a0a-9f81-89ae382aaf45
# ╟─2136ca01-87b4-47e9-93b1-25892896413b
# ╠═e664515e-d86b-4a91-9d47-2d77ad19d06a
# ╠═e2f12441-c8e7-486c-99b4-0ba302d41895
# ╟─5bcd4573-7656-4ebe-ba88-45adf2cc4127
# ╠═a25db436-2022-4c33-87fc-36f9b507823f
# ╠═b4433a1c-7d03-440f-af08-51d6d5f6d822
# ╟─f52a6b52-8202-4c3f-a7b4-4280ffc191d4
# ╠═4da06cb3-5086-4fe8-b728-a3f8e7e889c7
# ╠═4cb723a0-d3de-4a98-8d41-0e045a161d07
# ╟─252f5036-62f7-45ef-b7ca-5865aab4cda8
# ╠═a5c6b75a-be4b-4c82-ab87-2b72463b6303
# ╠═068905e6-9cde-4cc9-a042-a36b79c1a234
# ╟─8a3d4776-5968-450d-982f-00249caa82ae
# ╠═bbd279b4-184d-4e6f-bacf-b468d5c1fd36
# ╠═24a1c30c-a21c-471d-bcf7-d90fe7200a8c
# ╟─afee73f3-740d-4114-96c5-b50b1e025e6c
# ╠═d3996cd9-6890-4445-b2fe-c38ecc76ce77
# ╠═26d9f8ea-5d85-4cdd-b500-0e2d63f60e81
# ╠═13ce095f-2ba2-42a5-b753-cd6a43f024af
# ╠═7310762a-3c02-42e0-bbb7-85011ba4e144
