### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 0cd5e530-b1ae-11f1-b3de-a3d637e4cd5c
using Pkg; Pkg.activate("..")

# ╔═╡ d84fce92-1c9c-4ebd-8650-08e2087a8123
using InteractiveUtils

# ╔═╡ a2815bea-2f98-4067-b420-731c7c07dd71
using Plots

# ╔═╡ 6c6cbf81-1ee0-44bb-b247-b38057081b87
using ModelingToolkit

# ╔═╡ 77bbdf6a-67bf-4787-acae-42c6cc4ddb76
using ModelingToolkit: t_nounits as t, D_nounits as D

# ╔═╡ bc2755ba-0f59-4c63-ab42-81186500bd8f
using OrdinaryDiffEq

# ╔═╡ 06cbfa7b-c852-410a-a86f-57fec5f13068
using PlutoUI

# ╔═╡ 67a173da-d753-4444-9adc-4144cc8c2f1f
md"""
# Glass emptying
"""

# ╔═╡ 5f1b78ac-b41c-496e-9f99-0f7848690e88
md"""
Here is the setup of the problem:

$(LocalResource(joinpath(@__DIR__, "figures", "glass_filling_drawing.png"), :width => 400))
"""

# ╔═╡ ca001eb4-877e-49c0-894d-32ac9499d3b5
md"""
```math
\frac{dV}{dt} = A(h)\frac{dh}{dt} = -Q
\quad\Rightarrow\quad
\frac{dh}{dt} = -\frac{Q}{\pi (1.5 + 0.1h)^2}
```
"""

# ╔═╡ f56ab0d7-e618-44dc-9022-4ebcb267018d
@variables V(t) h(t)

# ╔═╡ f4f8f678-b198-4322-84db-695752dc7f46
@parameters Q = 0.1

# ╔═╡ 9d12ebc4-c0b8-4dc1-a99f-2628a4b13c2a
h0 = 5 # m

# ╔═╡ bdc248a8-559c-4c58-aca5-a8e930b5ce62
V0 = pi*((1.5+2)/2)^2*h0 # m³

# ╔═╡ 3bbb8572-bb3d-483d-8ef0-46e9b8e649ae
eqns = [D(V) ~ -Q,
	   D(h) ~ -Q/(pi*(1.5+0.1*h)^2)]

# ╔═╡ 316ae129-5418-482e-bdbb-51529e383426
@mtkbuild glass_sys = ODESystem(eqns, t)

# ╔═╡ 5875e4f3-4b90-4133-b5e9-8c1cd9f7e98a
glass_prob = ODEProblem(glass_sys, [V => V0, h =>h0], [0, 500])

# ╔═╡ 6c80922c-7533-471f-844c-9d50cffca030
glass_sol = solve(glass_prob, Tsit5());

# ╔═╡ 9d4d26f4-772b-4bce-a21a-918fa6a1bad2
plot(glass_sol, idxs = [h], label = "Height", xlabel ="Time (s)",
	ylabel = "height (m)", title = "Change of height")

# ╔═╡ Cell order:
# ╠═0cd5e530-b1ae-11f1-b3de-a3d637e4cd5c
# ╠═d84fce92-1c9c-4ebd-8650-08e2087a8123
# ╠═a2815bea-2f98-4067-b420-731c7c07dd71
# ╠═6c6cbf81-1ee0-44bb-b247-b38057081b87
# ╠═77bbdf6a-67bf-4787-acae-42c6cc4ddb76
# ╠═bc2755ba-0f59-4c63-ab42-81186500bd8f
# ╠═06cbfa7b-c852-410a-a86f-57fec5f13068
# ╟─67a173da-d753-4444-9adc-4144cc8c2f1f
# ╟─5f1b78ac-b41c-496e-9f99-0f7848690e88
# ╟─ca001eb4-877e-49c0-894d-32ac9499d3b5
# ╠═f56ab0d7-e618-44dc-9022-4ebcb267018d
# ╠═f4f8f678-b198-4322-84db-695752dc7f46
# ╠═9d12ebc4-c0b8-4dc1-a99f-2628a4b13c2a
# ╠═bdc248a8-559c-4c58-aca5-a8e930b5ce62
# ╠═3bbb8572-bb3d-483d-8ef0-46e9b8e649ae
# ╠═316ae129-5418-482e-bdbb-51529e383426
# ╠═5875e4f3-4b90-4133-b5e9-8c1cd9f7e98a
# ╠═6c80922c-7533-471f-844c-9d50cffca030
# ╠═9d4d26f4-772b-4bce-a21a-918fa6a1bad2
