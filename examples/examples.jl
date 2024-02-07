### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 4952d2aa-6b42-11ec-15ea-ed5764c1810b
begin
	import Pkg; Pkg.activate()
	using Revise
	using Controlz, DataFrames, CairoMakie, ColorSchemes, Printf
end

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

# ╔═╡ a6b0f813-91ed-460b-b65f-76fdbd799e39
md"## numerically invert the Laplace transform of a function in the frequency domain, back into the time domain.
(usually an input)

$\mathcal{L}[t \cos(at)]= \dfrac{s^2-a^2}{(s^2+a^2)^2}$
"

# ╔═╡ 8a510cf6-a0a1-490c-a2fa-f17251c098ce
function invert_lt()
	a = π
	U = (s^2 - a^2) / (s^2 + a^2) ^ 2
	data = simulate(U, 8.0, nb_time_points=300)

	fig = viz_response(data, title="inverting an input U*(s)", ylabel="u*(t)",
	                   savename="../docs/src/tcosat.png")
	return fig
end

# ╔═╡ b1adc640-477f-42c7-92fa-d7de44f5fb8e
invert_lt()

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
	g = 3 / (s + 1)

	fig = bode_plot(g, savename="../docs/src/example_bode.png")
	return fig
end

# ╔═╡ 26d9f8ea-5d85-4cdd-b500-0e2d63f60e81
draw_bode_plot()

# ╔═╡ cb5572cd-4c1d-4291-bfa3-4b83f9c9c784
md"## servo response of a control system"

# ╔═╡ 08722d92-f2cd-4ec1-bfca-e21ebdb47ad7
function sim_servo_response()
	Kc = 1.0
	τI = 1.0
	pic = PIController(1.0, 1.0) 
	gc = TransferFunction(pic) # controller transfer functio

	gp = 3 / (4 * s + 1) # process transfer function

	g_ol = gc * gp # open-loop transfer function

	g_servo = g_ol / (1 + g_ol) # transfer function for servo response

	Y_sp = 1 / s # unit step set point change

	Y = g_servo * Y_sp # resulting output

	E = Y_sp - Y # error signal

	U = gc * E # resulting controller output

	# break output into P-action, I-action
	U_Paction = Kc * E
	U_Iaction = Kc / (τI * s) * E
	
	# simulate for y, u, ysp in the time domain
	final_time = 12.0
	y_data = simulate(Y, final_time)
	u_data = simulate(U, final_time)
	ysp_data = simulate(Y_sp, final_time)
	u_Paction_data = simulate(U_Paction, final_time)
	u_Iaction_data = simulate(U_Iaction, final_time)
	
	fig = Figure(resolution=(800, 600))
	axs = [Axis(fig[1, 1], ylabel="system output y*(t)"),
		   Axis(fig[2, 1], xlabel="time, t", ylabel="controller output u*(t)")
	]
	
	lines!(axs[1], y_data[:, :t], y_data[:, :output], label="response y*(t)")
	lines!(axs[1], ysp_data[:, :t], ysp_data[:, :output], linestyle=:dash, label="set point yₛₚ(t)")
	axislegend(axs[1], position=:rb)
	
	lines!(axs[2], u_data[:, :t], u_data[:, :output], label="total")
	lines!(axs[2], u_Paction_data[:, :t], u_Paction_data[:, :output], linestyle=:dash, label="P-action", linewidth=2)
	lines!(axs[2], u_Iaction_data[:, :t], u_Iaction_data[:, :output], linestyle=:dot, label="I-action", linewidth=2)
	axislegend(axs[2])

	Controlz.draw_axes.(axs)
	save("../docs/src/simple_servo_response.png", fig)
	return fig
end

# ╔═╡ b87d02fc-2c8c-4536-b244-f7573e11ac4e
sim_servo_response()

# ╔═╡ 3f41f47d-68ef-4845-96ff-03f62f0a4652
cool_theme[:palette]

# ╔═╡ 80e11f0b-f06c-4abb-b1b7-45a47032b8b8
md"## closed loop servo response with time delay"

# ╔═╡ d9b1ad99-adcb-49e2-be52-2f194400dec9
function sim_servo_response_w_delay()
	
	# PI controller transfer function
	pic = PIController(1.0, 3.0)
	gc = TransferFunction(pic)
	
	# process, sensor dynamics
	gu = 2 / (4 * s + 1) * exp(-0.5 * s)
	gm = 1 / (s + 1) * exp(-0.75 * s)
	gd = 6 / (6 * s + 1)
	
	# open-loop transfer function
	g_ol = gc * gu * gm
	
	# closed-loop transfer function for regulator response
	gr = ClosedLoopTransferFunction(gd, g_ol)
	
	# closed-loop transfer function for servo response
	gs = ClosedLoopTransferFunction(gu * gc, g_ol)
	
	# response to unit set point change
	Y = gs / s
	data = simulate(Y, 50.0)
	
	# # response to unit step in disturbance d
	# Y = gr / s
	# data = simulate(Y, 45.0)
	
	fig = viz_response(data, title="closed-loop servo response", 
		               savename="../docs/src/closed_loop_servo_time_delay.png")
	return fig
end

# ╔═╡ 6a68842e-ecb7-45d3-86b1-2744752ad3b0
sim_servo_response_w_delay()

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

# ╔═╡ ed90dd6e-61cf-4736-a9c4-b8ed9702e7bb
begin
	local fig = Figure()
	local ax = Axis(fig[1, 1], xlabel="hi cory, y*(t)", xlabelfont="Space Mono")
	fig
end

# ╔═╡ Cell order:
# ╠═4952d2aa-6b42-11ec-15ea-ed5764c1810b
# ╠═ecab99c5-d929-4a0a-9f81-89ae382aaf45
# ╟─2136ca01-87b4-47e9-93b1-25892896413b
# ╠═e664515e-d86b-4a91-9d47-2d77ad19d06a
# ╠═e2f12441-c8e7-486c-99b4-0ba302d41895
# ╟─5bcd4573-7656-4ebe-ba88-45adf2cc4127
# ╠═a25db436-2022-4c33-87fc-36f9b507823f
# ╠═b4433a1c-7d03-440f-af08-51d6d5f6d822
# ╟─a6b0f813-91ed-460b-b65f-76fdbd799e39
# ╠═8a510cf6-a0a1-490c-a2fa-f17251c098ce
# ╠═b1adc640-477f-42c7-92fa-d7de44f5fb8e
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
# ╟─cb5572cd-4c1d-4291-bfa3-4b83f9c9c784
# ╠═08722d92-f2cd-4ec1-bfca-e21ebdb47ad7
# ╠═b87d02fc-2c8c-4536-b244-f7573e11ac4e
# ╠═3f41f47d-68ef-4845-96ff-03f62f0a4652
# ╟─80e11f0b-f06c-4abb-b1b7-45a47032b8b8
# ╠═d9b1ad99-adcb-49e2-be52-2f194400dec9
# ╠═6a68842e-ecb7-45d3-86b1-2744752ad3b0
# ╟─2de0d044-1bf6-489f-9e54-026100fd7655
# ╠═b6a45495-e081-4e02-8a0e-cfe1d3df36ba
# ╠═8acaad0e-951d-4bc9-add8-19315bb9f4e6
# ╠═ed90dd6e-61cf-4736-a9c4-b8ed9702e7bb
