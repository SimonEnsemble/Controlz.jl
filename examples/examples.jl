### A Pluto.jl notebook ###
# v0.17.4

using Markdown
using InteractiveUtils

# ╔═╡ 4952d2aa-6b42-11ec-15ea-ed5764c1810b
import Pkg; Pkg.activate()

# ╔═╡ 5da3e954-e940-45a7-b355-45e78b5ce38b
using Controlz, DataFrames, CairoMakie, ColorSchemes, Colors, Printf

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
fig = draw_root_locus_diagram()

# ╔═╡ c7739f1b-f3d9-4473-91fe-947f29b499f0
current_axis().limits

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

# ╔═╡ 2de0d044-1bf6-489f-9e54-026100fd7655
md"## draw a `.gif` of a response"

# ╔═╡ b6a45495-e081-4e02-8a0e-cfe1d3df36ba
function test_make_gif()
	K = 2.0 # gain
	τ = 4.0 # time constant
	θ = 1.5 # time delay
	
	g = K * exp(-θ * s) / (τ * s + 1) # FOPTD transfer function

	U = 1 / s # step input
	Y = g * U # response

	t_final = 15.0
	data = simulate(Y, t_final)
	mk_gif(data)
end

# ╔═╡ 8acaad0e-951d-4bc9-add8-19315bb9f4e6
test_make_gif()

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
# ╠═c7739f1b-f3d9-4473-91fe-947f29b499f0
# ╟─afee73f3-740d-4114-96c5-b50b1e025e6c
# ╠═d3996cd9-6890-4445-b2fe-c38ecc76ce77
# ╠═26d9f8ea-5d85-4cdd-b500-0e2d63f60e81
# ╟─2de0d044-1bf6-489f-9e54-026100fd7655
# ╠═b6a45495-e081-4e02-8a0e-cfe1d3df36ba
# ╠═8acaad0e-951d-4bc9-add8-19315bb9f4e6
