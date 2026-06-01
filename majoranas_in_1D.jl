### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ fef75640-5d9f-11f1-9745-0788a5745753
begin
	using LinearAlgebra
	using Plots
	using Plots.PlotMeasures          # unidades mm para los márgenes
	using PlutoUI
	using Printf
	gr()
	default(framestyle = :box, gridalpha = 0.25, legendfontsize = 8,
		titlefontsize = 10, dpi = 110,
		left_margin = 6mm, bottom_margin = 6mm,   # evita que se corten los labels
		right_margin = 4mm, top_margin = 3mm)
	md"Paquetes y estilo cargados ✓"
end

# ╔═╡ fef64cb4-5d9f-11f1-86d6-b9cf870cb680
md"""
# Modos de Majorana en la cadena de Kitaev

## El modelo de Kitaev

El modelo de Kitaev (Kitaev, 2001) es el ejemplo mínimo y más transparente de un superconductor topológico unidimensional, y por ello el punto de partida ideal para entender de dónde salen los famosos *modos de Majorana*. Describe una cadena de ``N`` sitios sobre la que se mueven fermiones sin espín (*spinless*) mediante tres ingredientes físicos: un término de salto (*hopping*) de amplitud ``t`` que permite a los electrones moverse entre sitios vecinos, un potencial químico ``\mu`` que fija cuántas partículas hay y el coste energético de ocupar cada sitio, y un emparejamiento superconductor de tipo ``p``-wave de amplitud ``\Delta`` que crea y destruye electrones por pares en sitios contiguos. Su Hamiltoniano es

```math
H = -\mu \sum_{j=1}^{N} c_j^\dagger c_j \;-\; \sum_{j=1}^{N-1}\left( t\, c_j^\dagger c_{j+1} + \Delta\, c_j c_{j+1} + \text{h.c.} \right),
```

donde ``c_j^\dagger`` y ``c_j`` crean y aniquilan un fermión en el sitio ``j``. La clave del modelo es reescribir cada fermión ordinario como una superposición de **dos operadores de Majorana** ``\gamma_{2j-1}`` y ``\gamma_{2j}`` (operadores hermíticos que son su propia antipartícula, ``\gamma = \gamma^\dagger``). En el régimen topológico, el Hamiltoniano empareja los Majoranas de sitios *vecinos*, dejando **desemparejados** un Majorana en cada extremo de la cadena: estos dos modos de borde, espacialmente separados, se combinan en un único estado fermiónico de energía cero, robusto frente al desorden y protegido por la topología del sistema.

A lo largo de este notebook construiremos esta física paso a paso: (1) veremos **cómo emergen los modos de Majorana** al reescribir el Hamiltoniano en la base de Majorana; (2) los **visualizaremos** representando la función de onda del modo de energía cero localizada en los extremos; (3) exploraremos los **distintos regímenes** que se obtienen al variar ``t``, ``\mu`` y ``\Delta``, distinguiendo la fase trivial (``|\mu| > 2|t|``) de la fase topológica (``|\mu| < 2|t|``, con ``\Delta \neq 0``); (4) calcularemos el **Pfaffiano** del Hamiltoniano en la base de Majorana como invariante topológico ``\mathbb{Z}_2`` que distingue ambas fases sin ambigüedad; y (5) conectaremos con el experimento estudiando el **pico de conductancia a voltaje cero** (*zero-bias conductance peak*), la firma experimental más característica de un modo de Majorana en una medida de transporte.
"""

# ╔═╡ fef755fa-5d9f-11f1-87ad-230064c776e7
md"""
## 0. Preparativos

Usamos `LinearAlgebra` (diagonalización), `Plots` (figuras) y `PlutoUI` (los **controles
interactivos** con *sliders*). Trabajaremos con la representación de Bogoliubov-de Gennes (BdG) de la
cadena de Kitaev y, más adelante, con su forma en la base de Majorana.

> **Nota:** este es un notebook **reactivo** de Pluto. Mueve los *sliders* y las figuras se
> recalculan al instante. Todo es **Julia puro**: incluso el cálculo de transporte de la
> sección 10 está implementado de forma nativa, así que no hay dependencias externas.
"""

# ╔═╡ fef7564a-5d9f-11f1-881a-4f94b82b21c4
TableOfContents(title = "📖 Contenido", depth = 2)

# ╔═╡ fef7567c-5d9f-11f1-9129-7730975337cb
md"""
## 1. El Hamiltoniano de Bogoliubov-de Gennes

Escribimos el Hamiltoniano en la base de Nambu
``\Psi = (c_1,\dots,c_N,\,c_1^\dagger,\dots,c_N^\dagger)^T`` como
``H = \tfrac12\,\Psi^\dagger\, H_{\rm BdG}\,\Psi + \text{cte}``, con

```math
H_{\rm BdG} = \begin{pmatrix} H_0 & D \\ -D^{*} & -H_0^{*}\end{pmatrix},
```

donde ``H_0`` es la parte normal (``-\mu`` en la diagonal, ``-t`` en los saltos a primeros
vecinos) y ``D`` es la matriz **antisimétrica** de emparejamiento (``D_{j,j+1}=\Delta``,
``D_{j+1,j}=-\Delta``). El factor ``\tfrac12`` y la duplicación de grados de libertad
(electrón + hueco) son el precio de tratar la superconductividad: el espectro resulta
simétrico ``E \leftrightarrow -E`` (simetría partícula-hueco).
"""

# ╔═╡ fef75690-5d9f-11f1-af3c-391647f08a24
"Matriz BdG 2N × 2N de la cadena de Kitaev abierta (cadena finita)."
function kitaev_bdg(N, t, μ, Δ)
	H0 = zeros(ComplexF64, N, N)
	D  = zeros(ComplexF64, N, N)
	for j in 1:N
		H0[j, j] = -μ                       # potencial químico
	end
	for j in 1:N-1
		H0[j, j+1] = -t;  H0[j+1, j] = -t   # hopping
		D[j, j+1]  =  Δ;  D[j+1, j]  = -Δ   # pairing (antisimétrico)
	end
	return [H0  D; -conj(D)  -conj(H0)]
end

# ╔═╡ fef75698-5d9f-11f1-8fae-a3c138b5a123
let
	# comprobación rápida: la matriz BdG es hermítica
	H = kitaev_bdg(4, 1.0, 0.5, 0.8)
	"Hermítica: $(isapprox(H, H'))"
end

# ╔═╡ fef756a4-5d9f-11f1-961f-cb221fc79870
md"""
## 2. Espectro: la aparición del modo de energía cero

Diagonalizamos ``H_{\rm BdG}`` para una cadena en la **fase topológica** (``|\mu|<2t``) y en
la **fase trivial** (``|\mu|>2t``). En la fase topológica aparece un par de estados
*pegados a ``E=0``* (separados del resto por el gap): son los dos modos de Majorana de los
extremos combinados. En la fase trivial no hay nada en el gap.
"""

# ╔═╡ fef756ae-5d9f-11f1-bc09-f5d1ea91bbeb
let
	N, t, Δ = 25, 1.0, 1.0
	plots = []
	for (μ, title) in [(0.5, "Topológica  |μ|<2t"), (3.0, "Trivial  |μ|>2t")]
		E = eigvals(Hermitian(kitaev_bdg(N, t, μ, Δ)))
		p = scatter(E; ms = 4, legend = false, title = "$title  (μ=$μ)",
			xlabel = "índice del autovalor", ylabel = "Energía E")
		hline!(p, [0]; ls = :dash, c = :black, lw = 0.6)
		push!(plots, p)
	end
	plot(plots...; layout = (1, 2), size = (900, 380))
end

# ╔═╡ fef756b8-5d9f-11f1-a111-279fef9838e0
md"""
**Explorador interactivo del espectro:** mueve ``\mu`` y ``\Delta`` y observa cómo el modo
se pega a ``E=0`` dentro de la fase topológica.

``\mu`` = $(@bind μ_sp Slider(-4:0.1:4, default = 0.5, show_value = true))

``\Delta`` = $(@bind Δ_sp Slider(0:0.1:2, default = 1.0, show_value = true))

``N`` = $(@bind N_sp Slider(4:60, default = 25, show_value = true))
"""

# ╔═╡ fef756cc-5d9f-11f1-9a43-819a539846c1
let
	E = eigvals(Hermitian(kitaev_bdg(N_sp, 1.0, μ_sp, Δ_sp)))
	topo = abs(μ_sp) < 2.0 && Δ_sp != 0
	estado = topo ? "TOPOLÓGICA (modo a E=0)" : "trivial"
	scatter(E; ms = 4, c = (topo ? 2 : 1), legend = false, ylim = (-4, 4),
		xlabel = "índice del autovalor", ylabel = "Energía E", size = (680, 420),
		title = "μ=$(round(μ_sp,digits=2)), Δ=$(round(Δ_sp,digits=2))  →  $estado")
	hline!([0]; ls = :dash, c = :black, lw = 0.6)
end

# ╔═╡ fef756d6-5d9f-11f1-aaaf-778e2741c1a3
md"""
## 3. De fermiones a Majoranas

Cada fermión se descompone en dos operadores de Majorana hermíticos:
```math
\gamma_{2j-1} = c_j + c_j^\dagger, \qquad \gamma_{2j} = -i\,(c_j - c_j^\dagger),
\qquad \gamma_a = \gamma_a^\dagger,\ \{\gamma_a,\gamma_b\}=2\delta_{ab}.
```
En esta base el Hamiltoniano toma la forma ``H=\tfrac{i}{4}\sum_{ab}A_{ab}\,\gamma_a\gamma_b``
con ``A`` **real y antisimétrica** (``A=-A^{T}``). La transformación es
``A = \tfrac{i}{4}\,\Omega\,H_{\rm BdG}\,\Omega^\dagger`` (salvo factores de convención),
con ``\Omega`` la matriz unitaria que pasa de ``(c,c^\dagger)`` a ``(\gamma_{\rm impar},\gamma_{\rm par})``.

En el **punto especialmente didáctico** ``\mu=0,\ t=\Delta``, el Hamiltoniano empareja
``\gamma_{2j}`` con ``\gamma_{2j+1}`` (Majoranas de sitios *vecinos*) y deja completamente
libres ``\gamma_1`` (extremo izquierdo) y ``\gamma_{2N}`` (extremo derecho): dos modos de
energía exactamente cero, deslocalizados entre los dos bordes.
"""

# ╔═╡ fef756e0-5d9f-11f1-8626-d5fcfc82545e
"Devuelve la matriz real antisimétrica A en la base de Majorana."
function to_majorana(Hbdg)
	twoN = size(Hbdg, 1)
	N = twoN ÷ 2
	# Ω: (c, c†) -> (γ_impar, γ_par)
	I_N = Matrix{ComplexF64}(I, N, N)
	Ω = zeros(ComplexF64, twoN, twoN)
	Ω[1:2:twoN, 1:N]     .=  I_N            # γ_{2j-1} = c_j + c_j†
	Ω[1:2:twoN, N+1:end] .=  I_N
	Ω[2:2:twoN, 1:N]     .=  im * I_N       # γ_{2j}   = -i(c_j - c_j†)
	Ω[2:2:twoN, N+1:end] .= -im * I_N
	Ω ./= sqrt(2)
	A = 0.5im * (Ω * Hbdg * Ω')
	Ar = real(A)
	return 0.5 * (Ar - transpose(Ar))       # antisimetriza por seguridad numérica
end

# ╔═╡ fef756ea-5d9f-11f1-99ae-6921d3e2e4f9
let
	A = to_majorana(kitaev_bdg(4, 1.0, 0.0, 1.0))
	nz = [i for i in 1:size(A, 1) if all(abs.(A[i, :]) .< 1e-9)]
	(
		es_real_antisimetrica = isapprox(A, -transpose(A)),
		filas_columnas_nulas = nz,
		interpretacion = "dos Majoranas desemparejados, uno en cada extremo (E=0)",
		A = round.(A; digits = 3),
	)
end

# ╔═╡ fef7570a-5d9f-11f1-b08c-0b41129fc141
md"""
## 4. Visualización de los modos de Majorana

Tomamos los dos autoestados de BdG con ``|E|`` más pequeño y construimos las dos
combinaciones de Majorana ``\gamma_\pm``. Representamos ``|\gamma(j)|^2`` sitio a sitio:
en la fase topológica cada modo está **localizado exponencialmente en un extremo**, y la
longitud de localización crece al acercarse a la transición.
"""

# ╔═╡ fef75712-5d9f-11f1-880f-615c8eb0e4fc
"|peso|² por sitio de los dos modos de Majorana de los extremos."
function majorana_wavefunctions(N, t, μ, Δ)
	F = eigen(Hermitian(kitaev_bdg(N, t, μ, Δ)))
	w, v = F.values, F.vectors
	# subespacio 2N×2 generado por los dos estados más cercanos a E=0
	Q = v[:, sortperm(abs.(w))[1:2]]
	# Los dos autoestados de BdG están deslocalizados sobre AMBOS extremos (son la
	# combinación fermiónica de los dos Majoranas) y el solver los devuelve con fases
	# arbitrarias; en el punto t=Δ, μ=0 son además exactamente degenerados. Para extraer
	# los Majoranas (uno por extremo) diagonalizamos el operador posición restringido al
	# subespacio: sus autovectores son las combinaciones izquierda / derecha.
	sites = vcat(1:N, 1:N)                  # sitio de cada fila (electrón y hueco)
	X = Hermitian(Q' * (sites .* Q))        # posición proyectada (2×2)
	c = eigen(X).vectors
	γL, γR = Q * c[:, 1], Q * c[:, 2]       # localizados en extremos opuestos
	dens(g) = abs2.(g[1:N]) .+ abs2.(g[N+1:end])
	return dens(γL), dens(γR)
end

# ╔═╡ fef75726-5d9f-11f1-9dbc-0be9347a1a2b
let
	N = 40
	plots = []
	for (μ, title) in [(0.0, "Profundo en la fase topológica (μ=0)"),
		(1.8, "Cerca de la transición (μ=1.8)")]
		γL, γR = majorana_wavefunctions(N, 1.0, μ, 1.0)
		sites = 1:N
		p = plot(sites, γL; m = :circle, ms = 3, label = "γ_izq",
			xlabel = "sitio j", ylabel = "|γ(j)|²", title = title)
		plot!(p, sites, γR; m = :square, ms = 3, label = "γ_der")
		push!(plots, p)
	end
	plot(plots...; layout = (1, 2), size = (1000, 380))
end

# ╔═╡ fef75730-5d9f-11f1-bd62-3de0a8163c23
md"""
**Explorador interactivo de la localización en los extremos:**

``\mu`` = $(@bind μ_ed Slider(-4:0.1:4, default = 0.0, show_value = true))

``\Delta`` = $(@bind Δ_ed Slider(0:0.1:2, default = 1.0, show_value = true))

``N`` = $(@bind N_ed Slider(4:80, default = 40, show_value = true))
"""

# ╔═╡ fef7573c-5d9f-11f1-ad24-93beaebf9fb0
let
	γL, γR = majorana_wavefunctions(N_ed, 1.0, μ_ed, Δ_ed)
	sites = 1:N_ed
	topo = abs(μ_ed) < 2.0 && Δ_ed != 0
	estado = topo ? "TOPOLÓGICA: un modo en cada borde" : "trivial: deslocalizado"
	plot(sites, γL; m = :circle, ms = 3, label = "γ_izq", size = (720, 420),
		xlabel = "sitio j", ylabel = "|γ(j)|²",
		title = "μ=$(round(μ_ed,digits=2))  →  $estado")
	plot!(sites, γR; m = :square, ms = 3, label = "γ_der")
end

# ╔═╡ fef75744-5d9f-11f1-b376-f32645221222
md"""
## 5. Función espectral y densidad local de estados (LDOS)

Para ver *a la vez* **dónde** y **a qué energía** viven los estados usamos la **función
espectral local de BdG**, que suma las componentes de electrón ``u_n`` y de hueco ``v_n``:
```math
A(j, E) = \sum_n \big(|u_n(j)|^2 + |v_n(j)|^2\big)\, \frac{1}{\pi}\frac{\eta}{(E-E_n)^2+\eta^2},
```
con ``\eta`` un ensanchamiento. La integral en energía da la LDOS del sitio ``j``.

> **Simetría partícula-hueco.** Cada autoestado ``(u_n,v_n)`` a energía ``E_n`` tiene una
> pareja ``(v_n^*,u_n^*)`` a ``-E_n``. Por eso la función espectral *electrónica*
> ``A_e=\sum_n|u_n|^2 L`` **no** es simétrica por sí sola: cumple ``A_e(j,E)=A_h(j,-E)``
> (el peso de electrón a ``+E`` es el de hueco a ``-E``). Es la suma ``A=A_e+A_h`` la que
> es **simétrica**, ``A(j,E)=A(j,-E)``, igual que el espectro de la sec. 2. Aquí dibujamos
> esa ``A`` simétrica; en la sec. 8 usaremos el canal electrónico ``A_e`` solo, porque es lo
> que mide un contacto de túnel de electrones.

En el mapa ``(E, j)``: en la **fase topológica** aparece una **línea brillante a ``E=0``
concentrada en los dos extremos** — la huella del modo de Majorana, separado del continuo
por el gap. En la **fase trivial** no hay peso a ``E=0`` y el extremo está dentro del gap.
"""

# ╔═╡ fef7574e-5d9f-11f1-a49e-63a88d4af75e
"""
Función espectral local A(j, E): filas = sitios, columnas = energía.

`canal = :total` usa el peso completo de BdG ``|u_n(j)|^2 + |v_n(j)|^2`` (electrón + hueco):
es **simétrico partícula-hueco**, ``A(j,E)=A(j,-E)``, como el espectro.
`canal = :electron` usa solo la componente electrón ``|u_n(j)|^2`` (lo que acopla un túnel de
electrones, sec. 8); ése **no** es simétrico, sino que cumple ``A_e(j,E)=A_h(j,-E)``.
"""
function ldos_map(N, t, μ, Δ, energies; η = 0.05, canal = :total)
	F = eigen(Hermitian(kitaev_bdg(N, t, μ, Δ)))
	w, v = F.values, F.vectors
	weight = canal === :electron ? abs2.(v[1:N, :]) :        # |u_n(j)|²
		abs2.(v[1:N, :]) .+ abs2.(v[N+1:end, :])             # |u_n|²+|v_n|² (N, 2N)
	lor = (η / π) ./ ((energies .- w').^2 .+ η^2)            # (E, 2N)
	return weight * lor'                                     # (N sitios, E)
end

# ╔═╡ fef75762-5d9f-11f1-a639-1b0db8fa8e2d
let
	N = 60
	E = range(-3, 3, length = 400)
	plots = []
	for (μ, title) in [(0.5, "Topológica (μ=0.5)"), (3.0, "Trivial (μ=3.0)")]
		A = ldos_map(N, 1.0, μ, 1.0, E; η = 0.06)
		p = heatmap(E, 1:N, A; c = :inferno, xlabel = "Energía E",
			ylabel = "sitio j", title = title, colorbar_title = "A(j,E)")
		push!(plots, p)
	end
	plot(plots...; layout = (1, 2), size = (1000, 400),
		plot_title = "Función espectral (simétrica E↔−E): el Majorana es la línea a E=0 en los bordes",
		plot_titlefontsize = 10)
end

# ╔═╡ fef75776-5d9f-11f1-91f8-ed0c37744c5d
md"""
**Explorador interactivo de la función espectral / LDOS:**

``\mu`` = $(@bind μ_ld Slider(-4:0.1:4, default = 0.5, show_value = true))

``\Delta`` = $(@bind Δ_ld Slider(0:0.1:2, default = 1.0, show_value = true))

``\eta`` = $(@bind η_ld Slider(0.02:0.01:0.25, default = 0.06, show_value = true))

``N`` = $(@bind N_ld Slider(20:2:100, default = 60, show_value = true))
"""

# ╔═╡ fef75780-5d9f-11f1-927f-e3cc3b6416a1
let
	E = range(-3, 3, length = 300)
	A = ldos_map(N_ld, 1.0, μ_ld, Δ_ld, E; η = η_ld)
	topo = abs(μ_ld) < 2.0 && Δ_ld != 0
	heatmap(E, 1:N_ld, A; c = :inferno, size = (720, 460),
		xlabel = "Energía E", ylabel = "sitio j", colorbar_title = "A(j,E)",
		title = "μ=$(round(μ_ld,digits=2)), Δ=$(round(Δ_ld,digits=2))  →  $(topo ? "TOPOLÓGICA" : "trivial")")
end

# ╔═╡ fef75794-5d9f-11f1-8289-89f9d7298daa
md"""
## 6. Regímenes de ``t``, ``\mu`` y ``\Delta``

Dos barridos para entender el espacio de parámetros:

1. **Energía del modo de borde frente a ``\mu``**: la energía del estado más bajo se
   pega a cero dentro de ``|\mu|<2t`` y se despega fuera. La transición está en ``|\mu|=2t``,
   donde el *gap del bulk* se cierra.
2. **Diagrama de fases en el plano ``(\mu,\Delta)``**: el gap se cierra en las líneas
   ``|\mu|=2t`` (para ``\Delta\neq0``) y en ``\Delta=0`` (sin emparejamiento no hay fase
   topológica). El interior ``|\mu|<2t,\ \Delta\neq0`` es la región topológica.
"""

# ╔═╡ fef7579c-5d9f-11f1-9e8e-2ffee7f512e5
let
	N, t, Δ = 60, 1.0, 1.0
	mus = range(-4, 4, length = 241)
	E0  = Float64[]   # energía del estado más bajo (>0)
	gap = Float64[]   # gap del bulk (tercer estado positivo, aprox)
	for μ in mus
		w = sort(abs.(eigvals(Hermitian(kitaev_bdg(N, t, μ, Δ)))))
		push!(E0, w[1])
		push!(gap, w[3])
	end
	plot(mus, E0; label = "modo de borde |E₀|", lw = 2, size = (720, 420),
		xlabel = "μ / t", ylabel = "Energía",
		title = "Modo de borde y gap frente a μ  (t=Δ=1)")
	plot!(mus, gap; label = "gap del bulk", alpha = 0.7)
	vspan!([-2, 2]; c = 2, alpha = 0.10, label = "región topológica")
	vline!([-2, 2]; c = :red, ls = :dash, lw = 0.8, label = "")
	annotate!(0, 0.85 * maximum(gap), text("TOPOLÓGICA", :green, 9))
end

# ╔═╡ fef757a8-5d9f-11f1-a24c-fdc49fe87dfa
let
	N, t = 30, 1.0
	mus = range(-4, 4, length = 121)
	Δs  = range(-2, 2, length = 121)
	GAP = [minimum(abs.(eigvals(Hermitian(kitaev_bdg(N, t, μ, Δ))))) for Δ in Δs, μ in mus]
	heatmap(mus, Δs, GAP; c = :viridis, size = (680, 520), colorbar_title = "min|E|",
		xlabel = "μ / t", ylabel = "Δ / t",
		title = "Diagrama de fases: zonas oscuras = modo a energía cero")
	vline!([-2, 2]; c = :white, ls = :dash, lw = 0.8, label = "")
	hline!([0]; c = :white, ls = :dot, lw = 0.8, label = "")
end

# ╔═╡ fef757bc-5d9f-11f1-b730-3b4c20c4c8e9
md"""
## 7. Invariantes del bulk: el Pfaffiano, el winding y la fase de Zak

¿Cómo distinguir las dos fases *sin* mirar los bordes? Con un **invariante del bulk**.
Kitaev demostró que el número de Majorana
```math
\mathcal{M} = \operatorname{sgn}\big[\operatorname{Pf}B(0)\big]\,
              \operatorname{sgn}\big[\operatorname{Pf}B(\pi)\big] = \pm 1
```
clasifica la cadena: ``\mathcal{M}=-1`` es topológica (un Majorana por extremo) y
``\mathcal{M}=+1`` es trivial. Aquí ``B(k)`` es el Hamiltoniano de Bloch en la base de
Majorana en los dos momentos invariantes bajo partícula-hueco, ``k=0`` y ``k=\pi``.

El **Pfaffiano** de una matriz antisimétrica satisface ``\operatorname{Pf}(A)^2=\det(A)`` y
es, en esencia, una "raíz cuadrada del determinante" con signo bien definido.
Implementamos el algoritmo de Parlett-Reid (no basta con ``\sqrt{\det}``, que pierde el signo).
"""

# ╔═╡ fef757ce-5d9f-11f1-89d3-b5623a62ff4e
"Pfaffiano de una matriz antisimétrica real (algoritmo de Parlett-Reid)."
function pfaffian(Ain)
	A = Matrix{Float64}(Ain)
	n = size(A, 1)
	@assert size(A, 1) == size(A, 2) && iseven(n)
	pf = 1.0
	for k in 1:2:n-1
		# pivoteo: mayor |A[k+1:n, k]|
		q = k + argmax(abs.(A[k+1:n, k]))
		if q != k + 1
			A[[k+1, q], :] = A[[q, k+1], :]
			A[:, [k+1, q]] = A[:, [q, k+1]]
			pf = -pf
		end
		A[k+1, k] == 0 && return 0.0
		pf *= A[k, k+1]
		if k + 2 <= n
			τ = A[k, k+2:n] ./ A[k, k+1]
			col = A[k+2:n, k+1]                      # <- columna (no fila)
			A[k+2:n, k+2:n] .+= τ * col' - col * τ'
		end
	end
	return pf
end

# ╔═╡ fef757da-5d9f-11f1-aa5c-279ebcdd70e7
let
	# test: Pf² == det
	M = randn(6, 6); M = M - M'
	(Pf2 = pfaffian(M)^2, det = det(M))
end

# ╔═╡ fef757e4-5d9f-11f1-851d-a54019eca6b4
begin
	"B(k): Hamiltoniano de Bloch 2×2 de Kitaev en la base de Majorana (real antisim.)."
	function bloch_majorana(k, t, μ, Δ)
		# H(k) en base Nambu (c_k, c_{-k}†):  ξ(k) τ_z + 2Δ sin k τ_y
		ξ = -2t * cos(k) - μ
		Hk = [ξ              -2im * Δ * sin(k);
			  2im * Δ * sin(k)  -ξ]
		U = [1 1; im -im] ./ sqrt(2)              # Nambu -> Majorana
		return real(im * (U * Hk * U'))
	end

	"Número de Majorana M = ±1 (−1: topológica, +1: trivial)."
	majorana_number(t, μ, Δ) =
		sign(pfaffian(bloch_majorana(0.0, t, μ, Δ))) *
		sign(pfaffian(bloch_majorana(π, t, μ, Δ)))
end

# ╔═╡ fef757ee-5d9f-11f1-81bf-7501dcb720ad
let
	for μ in (0.0, 1.0, 3.0)
		@printf "μ=%.1f:  M = %+d\n" μ Int(majorana_number(1.0, μ, 1.0))
	end
	mus = range(-4, 4, length = 401)
	M = [majorana_number(1.0, μ, 1.0) for μ in mus]
	plot(mus, M; lw = 2, legend = false, size = (720, 330), left_margin = 14mm,
		xlabel = "μ / t", ylabel = "M",
		yticks = ([-1, 1], ["-1 (topológica)", "+1 (trivial)"]),
		title = "Invariante de Pfaffiano frente a μ  (t=Δ=1)")
	vline!([-2, 2]; c = :red, ls = :dash, lw = 0.8)
end

# ╔═╡ 45be7096-5da4-11f1-9944-1913e66db732
md"""
### Otro invariante del bulk: el número de *winding* y la fase de Zak

El Pfaffiano da un ``\mathbb{Z}_2`` (``\pm1``) evaluado en los puntos especiales ``k=0,\pi``.
Hay un invariante hermano, geométrico y muy visual, que usa **toda** la zona de Brillouin.

Con parámetros reales, la cadena de Kitaev tiene **simetría quiral** (``\{H(k),\tau_x\}=0``),
así que pertenece a la clase **BDI**, con un invariante entero ``\mathbb{Z}``: el **número
de winding** ``w``. Escribiendo el Hamiltoniano de Bloch en base de Nambu como
``H(k) = h_z(k)\,\tau_z + h_y(k)\,\tau_y`` con
```math
h_z(k) = -2t\cos k - \mu, \qquad h_y(k) = 2\Delta\sin k,
```
el vector ``(h_z,h_y)`` recorre una **elipse** centrada en ``(-\mu,0)`` con semiejes
``(2t,2\Delta)`` al variar ``k\in[0,2\pi)``. El número de winding cuenta **cuántas veces
rodea esa curva al origen**:
```math
w = \frac{1}{2\pi}\oint dk\,\frac{d}{dk}\arg\big[h_z(k)+i\,h_y(k)\big]
  = \begin{cases} \pm1 & |\mu|<2t\ (\text{topológica}) \\ 0 & |\mu|>2t\ (\text{trivial}).\end{cases}
```
``|w|`` es exactamente el **número de Majoranas por borde** (correspondencia *bulk-boundary*).

La **fase de Zak** (la holonomía de Berry sobre la BZ, es decir, el *Wilson loop* 1D) es la
versión "integral de línea" de lo mismo: la calculamos numéricamente discretizando la BZ,
```math
\gamma_{\rm Zak} = -\,\text{Im}\,\ln\prod_{j} \langle u(k_j)\,|\,u(k_{j+1})\rangle ,
```
con ``|u(k)\rangle`` el autoestado de la banda ocupada. Para este sistema quiral vale
``\gamma_{\rm Zak}=\pi`` en la fase topológica y ``0`` en la trivial (relación
``\gamma_{\rm Zak}=\pi\,w \bmod 2\pi``).

> **Salvedad D vs BDI.** El ``\mathbb{Z}`` del winding existe porque hay simetría quiral. Si
> añades un término que la rompe (p. ej. un *hopping* complejo genérico), el invariante
> colapsa al ``\mathbb{Z}_2`` del Pfaffiano: solo sobrevive la paridad de ``w``.
"""

# ╔═╡ 45bf77fc-5da4-11f1-9426-33509b8395ba
begin
	# H(k) = h_z τ_z + h_y τ_y  (Bloch BdG en base de Nambu)
	h_z(k, t, μ) = -2t * cos(k) - μ
	h_y(k, Δ)    = 2Δ * sin(k)
	bloch_nambu(k, t, μ, Δ) = [h_z(k,t,μ)        -im*h_y(k,Δ);
		                        im*h_y(k,Δ)   -h_z(k,t,μ)]

	"Número de winding w (clase BDI): vueltas de q(k)=h_z+i h_y alrededor del origen."
	function winding_number(t, μ, Δ; M = 2000)
		ks = range(0, 2π, length = M + 1)[1:end-1]
		q = [h_z(k,t,μ) + im * h_y(k,Δ) for k in ks]
		sum(angle(q[mod1(j+1, M)] / q[j]) for j in 1:M) / (2π)
	end

	"Fase de Zak (Wilson loop discreto sobre la banda ocupada), en [0, 2π)."
	function zak_phase(t, μ, Δ; M = 2000)
		ks = range(0, 2π, length = M + 1)[1:end-1]
		us = [eigen(Hermitian(bloch_nambu(k,t,μ,Δ))).vectors[:, 1] for k in ks]
		W = prod(dot(us[i], us[mod1(i+1, M)]) for i in 1:M)
		mod(-angle(W), 2π)
	end
end

# ╔═╡ 45bf784c-5da4-11f1-8a05-7de32ed803cb
let
	ks = range(0, 2π, length = 400)
	plots = []
	for (μ, title) in [(0.5, "Topológica (μ=0.5)"), (3.0, "Trivial (μ=3.0)")]
		hzs = [h_z(k, 1.0, μ) for k in ks]
		hys = [h_y(k, 1.0) for k in ks]
		w = round(Int, winding_number(1.0, μ, 1.0))
		z = zak_phase(1.0, μ, 1.0)
		p = plot(hzs, hys; lw = 2, legend = false, aspect_ratio = :equal,
			xlabel = "h_z(k)", ylabel = "h_y(k)",
			title = "$title\nw=$w,  Zak=$(round(z/π,digits=2))π")
		scatter!(p, [0], [0]; m = :star5, ms = 9, c = :red)          # el origen
		push!(plots, p)
	end
	plot(plots...; layout = (1, 2), size = (900, 430),
		plot_title = "¿La elipse (h_z, h_y) rodea el origen (★)? → w=1 topológica, w=0 trivial",
		plot_titlefontsize = 9)
end

# ╔═╡ 45bf786a-5da4-11f1-8ca0-09dc71bb9675
md"""
**Explorador interactivo del winding:** mueve ``\mu`` y mira cómo la elipse deja de rodear
al origen (★) justo en ``|\mu|=2t``.

``\mu`` = $(@bind μ_wd Slider(-4:0.1:4, default = 0.5, show_value = true))

``\Delta`` = $(@bind Δ_wd Slider(0.1:0.1:2, default = 1.0, show_value = true))
"""

# ╔═╡ 45bf7892-5da4-11f1-ab20-1f27e393a807
let
	ks = range(0, 2π, length = 400)
	hzs = [h_z(k, 1.0, μ_wd) for k in ks]
	hys = [h_y(k, Δ_wd) for k in ks]
	w = round(Int, winding_number(1.0, μ_wd, Δ_wd))
	z = zak_phase(1.0, μ_wd, Δ_wd)
	topo = abs(μ_wd) < 2.0
	plot(hzs, hys; lw = 2, legend = false, aspect_ratio = :equal, size = (560, 520),
		c = (topo ? 2 : 1), xlabel = "h_z(k)", ylabel = "h_y(k)",
		title = "μ=$(round(μ_wd,digits=2)), Δ=$(round(Δ_wd,digits=2))  →  " *
			"w=$w, Zak=$(round(z/π,digits=2))π  ($(topo ? "TOPOLÓGICA" : "trivial"))")
	scatter!([0], [0]; m = :star5, ms = 10, c = :red)
end

# ╔═╡ fef757f8-5d9f-11f1-b13f-dbc40b051232
md"""
## 8. El pico de conductancia a voltaje cero (ZBCP)

La firma experimental de un Majorana en transporte es un **pico de conductancia
diferencial a voltaje cero**. En el límite túnel, la conductancia que mide un contacto
acoplado al extremo de la cadena es proporcional a la **densidad local de estados (LDOS)
electrónica del extremo**:
```math
\frac{dI}{dV}(V) \;\propto\; A_{e}(\text{extremo},\,E=eV)
= \sum_n |u_n(\text{extremo})|^2\,\frac{1}{\pi}\frac{\eta}{(eV-E_n)^2+\eta^2}.
```
El modo de Majorana de energía cero produce un **pico centrado en ``V=0``** en la fase
topológica; en la fase trivial el extremo está dentro del gap y la conductancia se anula
a bajo voltaje. (Este es un proxy espectral honesto; un cálculo de transporte completo
—sección 10— da además la cuantización ``2e^2/h`` del pico.)
"""

# ╔═╡ fef75800-5d9f-11f1-8e03-4f7a65ada5ea
"dI/dV (límite túnel) ~ LDOS *electrónica* del extremo (canal de electrones, sec. 5)."
function end_ldos(N, t, μ, Δ, energies; η = 0.04, lado = :izq)
	A = ldos_map(N, t, μ, Δ, energies; η = η, canal = :electron)
	return lado === :izq ? A[1, :] : A[end, :]
end

# ╔═╡ fef75816-5d9f-11f1-95b6-df45280bd3b9
let
	V = range(-3, 3, length = 601)
	p = plot(size = (720, 420), xlabel = "voltaje  eV / t",
		ylabel = "dI/dV  (u.a., ∝ LDOS extremo)",
		title = "Pico de conductancia a voltaje cero (ZBCP) — proxy espectral")
	for (μ, lbl, c) in [(0.5, "Topológica (μ=0.5)", 1), (3.0, "Trivial (μ=3.0)", 2)]
		g = end_ldos(60, 1.0, μ, 1.0, V; η = 0.05)
		plot!(p, V, g; label = lbl, c = c, lw = 2)
	end
	vline!(p, [0]; c = :black, ls = :dot, lw = 0.6, label = "")
	p
end

# ╔═╡ fef75820-5d9f-11f1-95c3-91c9618303a1
md"""
**Explorador interactivo del ZBCP** (proxy LDOS del extremo):

``\mu`` = $(@bind μ_zb Slider(-4:0.1:4, default = 0.5, show_value = true))

``\Delta`` = $(@bind Δ_zb Slider(0:0.1:2, default = 1.0, show_value = true))

``N`` = $(@bind N_zb Slider(20:2:100, default = 60, show_value = true))
"""

# ╔═╡ fef7582a-5d9f-11f1-83a7-41107188065b
let
	V = range(-3, 3, length = 400)
	g = end_ldos(N_zb, 1.0, μ_zb, Δ_zb, V; η = 0.05)
	topo = abs(μ_zb) < 2.0 && Δ_zb != 0
	plot(V, g; c = (topo ? 1 : 2), lw = 2, legend = false, size = (720, 420),
		xlabel = "voltaje eV/t", ylabel = "dI/dV (u.a.)",
		title = "μ=$(round(μ_zb,digits=2))  →  $(topo ? "ZBCP presente (topológica)" : "sin pico (trivial)")")
	vline!([0]; c = :black, ls = :dot, lw = 0.6)
end

# ╔═╡ fef75834-5d9f-11f1-a38d-41975b7dcd99
md"""
## 9. Resumen

| Concepto | Diagnóstico en el notebook |
|---|---|
| Modos de Majorana | par de estados a ``E\approx0`` separados del gap (sec. 2) |
| Localización | ``|\gamma(j)|^2`` en cada extremo (sec. 4) |
| Función espectral / LDOS | línea a ``E=0`` en los bordes (sec. 5) |
| Regímenes ``t,\mu,\Delta`` | fase topológica ``|\mu|<2t,\ \Delta\neq0`` (sec. 6) |
| Invariante de bulk (``\mathbb{Z}_2``) | Pfaffiano ``\mathcal{M}=-1`` en la fase topológica (sec. 7) |
| Invariante de bulk (``\mathbb{Z}``) | número de winding ``w=1``, fase de Zak ``=\pi`` (sec. 7) |
| Firma experimental | pico de conductancia a ``V=0`` (sec. 8) |
| Transporte / Andreev | ``R_{he}=1`` y ZBCP cuantizado ``2e^2/h`` (sec. 10) |

**Ideas para extender la charla:** desorden y robustez del modo cero, *fermion parity* y
cómputo cuántico topológico (trenzado de Majoranas), y el modelo de nanohilo realista
(Lutchyn-Oreg: acoplo espín-órbita + campo Zeeman + superconductor ``s``-wave).
"""

# ╔═╡ fef7583e-5d9f-11f1-b70b-6199465e6f42
md"""
## 10. La conductancia cuantizada (transporte nativo en Julia)

El proxy de la sección 8 reproduce la *forma* del ZBCP, pero no su **valor cuantizado**.
Con un cálculo de transporte real acoplando un **contacto normal metálico** al extremo de
la cadena se obtiene la conductancia de Andreev
```math
G = \frac{e^2}{h}\big(N_e - R_{ee} + R_{he}\big),
```
y un único modo de Majorana fuerza el **valor universal** ``G(V\!=\!0) = 2e^2/h`` en la
fase topológica, frente a ``G\approx0`` en la trivial.

> **Transporte en Julia puro.** El cálculo de transporte está implementado de forma nativa,
> sin librerías externas: construimos el dispositivo BdG en base Nambu por sitio, acoplamos
> un *lead* normal 1D mediante su **autoenergía** (función de Green de superficie, que es
> analítica para una cadena 1D) y obtenemos la matriz de *scattering* con la relación de
> Fisher-Lee. Así no hay dependencias externas y el cálculo es transparente.
"""

# ╔═╡ fef75852-5d9f-11f1-a9af-b9ea690da007
begin
	"Dispositivo de Kitaev en base Nambu por sitio (2 orbitales: electrón/hueco)."
	function device_H(N, t, μ, Δ)
		σy = ComplexF64[0 -im; im 0]
		σz = ComplexF64[1 0; 0 -1]
		H = zeros(ComplexF64, 2N, 2N)
		hop = -t * σz + im * Δ * σy
		for j in 1:N
			H[2j-1:2j, 2j-1:2j] = -μ * σz
		end
		for j in 1:N-1
			H[2j-1:2j, 2j+1:2j+2] = hop
			H[2j+1:2j+2, 2j-1:2j] = hop'
		end
		return H
	end

	"Función de Green de superficie (escalar) de una cadena 1D semiinfinita (onsite 0, hopping τ)."
	function surface_g(E, τ; η = 1e-9)
		z = E + im * η
		disc = sqrt(Complex(z^2 - 4τ^2))         # τ² g² - z g + 1 = 0
		g1 = (z - disc) / (2τ^2)
		g2 = (z + disc) / (2τ^2)
		return imag(g1) <= 0 ? g1 : g2           # rama retardada: Im g <= 0
	end
end

# ╔═╡ fef7585c-5d9f-11f1-a52d-833c5e4f6150
"G/(e²/h) = Nₑ − R_ee + R_he para un contacto normal acoplado al extremo izquierdo."
function conductance(N, t, μ, Δ, E; t_lead = 1.0, t_c = 1.0, barrier = 2.0, η = 1e-9)
	σz = ComplexF64[1 0; 0 -1]
	H = device_H(N, t, μ, Δ)
	H[1:2, 1:2] += barrier * σz                  # barrera túnel en el sitio 1
	g = surface_g(E, t_lead; η = η)              # lead normal con μ=0 (bloques e/h iguales)
	V = -t_c * σz                                # acoplo dispositivo-lead
	Σ = V * (g * I(2)) * V'                       # autoenergía 2×2 en el sitio 1
	Σfull = zeros(ComplexF64, 2N, 2N)
	Σfull[1:2, 1:2] = Σ
	Gr = inv((E + im * η) * I - H - Σfull)        # función de Green retardada
	Γ = im * (Σ - Σ')
	Γh = sqrt(Hermitian(Γ))
	S = -I(2) + im * (Γh * Gr[1:2, 1:2] * Γh)     # Fisher-Lee (un solo lead)
	Ree = abs2(S[1, 1])                           # reflexión normal
	Rhe = abs2(S[2, 1])                           # reflexión de Andreev
	return real(1 - Ree + Rhe)                    # Nₑ = 1
end

# ╔═╡ fef75866-5d9f-11f1-88fe-d77f4212fec6
let
	# comprobación: ~2 e²/h en la fase topológica, ~0 en la trivial
	res = []
	for (μ, lbl) in [(0.5, "topológica"), (3.0, "trivial")]
		g = conductance(150, 1.0, μ, 1.0, 1e-4)
		push!(res, (fase = lbl, μ = μ, G_en_V0 = round(g; digits = 3)))
	end
	res
end

# ╔═╡ fef75872-5d9f-11f1-824f-13213e62ca60
let
	V = range(-1, 1, length = 241) .+ 1e-4       # desplazada para evitar E=0 exacto
	p = plot(size = (720, 420), xlabel = "voltaje eV/t", ylabel = "G  (e²/h)",
		title = "ZBCP cuantizado (transporte real, Julia puro)")
	for (μ, lbl, c) in [(0.5, "Topológica (μ=0.5)", 1), (3.0, "Trivial (μ=3.0)", 2)]
		g = [conductance(150, 1.0, μ, 1.0, e) for e in V]
		plot!(p, V, g; label = lbl, c = c, lw = 2)
	end
	hline!(p, [2]; c = :gray, ls = :dash, lw = 0.8, label = "2e²/h")
	vline!(p, [0]; c = :black, ls = :dot, lw = 0.6, label = "")
	p
end

# ╔═╡ 45d08718-5da4-11f1-8648-b95a4c5a96da
md"""
### ¿Qué estamos viendo? La reflexión de Andreev

El contacto normal de la izquierda inyecta **electrones**. A energías por debajo del gap
del superconductor no hay estados de cuasipartícula donde entrar, así que el electrón
**se refleja**, y solo hay dos posibilidades:

- **Reflexión normal** (electrón → electrón, probabilidad ``R_{ee}``): rebota tal cual; no
  transfiere carga neta al superconductor.
- **Reflexión de Andreev** (electrón → **hueco**, probabilidad ``R_{he}``): el electrón se
  empareja con otro y entra como un **par de Cooper**, devolviendo un hueco. Esto transfiere
  carga ``2e`` y por eso *aumenta* la conductancia:
  ```math
  G = \frac{e^2}{h}\big(N_e - R_{ee} + R_{he}\big), \qquad R_{ee}+R_{he}=1 .
  ```

La clave: **un modo de Majorana a ``E=0`` fuerza reflexión de Andreev perfecta y resonante**
a voltaje cero, ``R_{he}=1`` y ``R_{ee}=0``, de modo que
```math
G(V=0) = \frac{2e^2}{h}\quad(\text{cuantizado y universal}),
```
independientemente de lo opaca que sea la barrera del contacto — esa robustez es la firma.
En la fase trivial, a ``V=0`` el electrón simplemente rebota (``R_{ee}=1,\ R_{he}=0``) y
``G\to0``.

En la figura de abajo se ve directamente: ``R_{he}`` (Andreev) sube hasta **1** en ``V=0``
solo en la fase topológica, mientras que ``R_{ee}`` (normal) cae a 0. Eso *es* el pico
cuantizado a ``2e^2/h``.
"""

# ╔═╡ 45d087ae-5da4-11f1-b5d4-d75b6fb9092d
"Probabilidades de reflexión normal (R_ee) y de Andreev (R_he) en el extremo izquierdo."
function andreev_R(N, t, μ, Δ, E; t_lead = 1.0, t_c = 1.0, barrier = 2.0, η = 1e-9)
	σz = ComplexF64[1 0; 0 -1]
	H = device_H(N, t, μ, Δ);  H[1:2, 1:2] += barrier * σz
	g = surface_g(E, t_lead; η = η);  V = -t_c * σz
	Σ = V * (g * I(2)) * V'
	Σfull = zeros(ComplexF64, 2N, 2N);  Σfull[1:2, 1:2] = Σ
	Gr = inv((E + im * η) * I - H - Σfull)
	Γh = sqrt(Hermitian(im * (Σ - Σ')))
	S = -I(2) + im * (Γh * Gr[1:2, 1:2] * Γh)
	return (Ree = abs2(S[1, 1]), Rhe = abs2(S[2, 1]))
end

# ╔═╡ 45d087d8-5da4-11f1-92de-6926d279b7df
let
	V = range(-1, 1, length = 241) .+ 1e-4
	plots = []
	for (μ, title) in [(0.5, "Topológica (μ=0.5)"), (3.0, "Trivial (μ=3.0)")]
		R = [andreev_R(150, 1.0, μ, 1.0, e) for e in V]
		p = plot(V, getfield.(R, :Rhe); lw = 2, label = "R_he (Andreev)", c = 2,
			xlabel = "voltaje eV/t", ylabel = "probabilidad", ylim = (-0.05, 1.05),
			title = title)
		plot!(p, V, getfield.(R, :Ree); lw = 2, label = "R_ee (normal)", c = 1)
		push!(plots, p)
	end
	plot(plots...; layout = (1, 2), size = (950, 400),
		plot_title = "Reflexión de Andreev perfecta (R_he→1) en V=0 ⟺ ZBCP cuantizado",
		plot_titlefontsize = 9)
end

# ╔═╡ fef7587a-5d9f-11f1-927f-b1b22b46a632
md"""
**Explorador interactivo del ZBCP real:** barre ``\mu`` a través de la transición.

``\mu`` = $(@bind μ_tr Slider(-4:0.1:4, default = 0.5, show_value = true))

barrera = $(@bind barrier_tr Slider(0:0.25:5, default = 2.0, show_value = true))
"""

# ╔═╡ fef7588e-5d9f-11f1-893a-1f654dc044e2
let
	V = range(-1, 1, length = 161) .+ 1e-4
	g = [conductance(150, 1.0, μ_tr, 1.0, e; barrier = barrier_tr) for e in V]
	topo = abs(μ_tr) < 2.0
	plot(V, g; c = (topo ? 1 : 2), lw = 2, legend = false, size = (720, 420),
		ylim = (-0.1, 2.4), xlabel = "voltaje eV/t", ylabel = "G (e²/h)",
		title = "μ=$(round(μ_tr,digits=2))  →  $(topo ? "ZBCP ≈ 2e²/h" : "sin pico")")
	hline!([2]; c = :gray, ls = :dash, lw = 0.8)
	vline!([0]; c = :black, ls = :dot, lw = 0.6)
end

# ╔═╡ fef75898-5d9f-11f1-81bb-cd27d653ddff
let
	# Mapa 2D G(V, μ): el ZBCP existe SÓLO dentro de la fase topológica |μ|<2t.
	# Usamos un broadening finito η>0: suaviza las resonancias Fabry-Pérot de la
	# cadena finita (que producen moteado cerca de μ=±2, donde se cierra el gap)
	# sin destruir la cuantización del ZBCP (plateau ≈ 1.93 e²/h).
	mus = range(-3.5, 3.5, length = 81)
	V = range(-1, 1, length = 81) .+ 1e-4
	Gmap = [conductance(120, 1.0, μ, 1.0, e; barrier = 2.0, η = 0.02) for e in V, μ in mus]
	heatmap(mus, V, Gmap; c = :inferno, clims = (0, 2), size = (740, 500),
		colorbar_title = "G (e²/h)", xlabel = "μ / t", ylabel = "voltaje eV/t",
		title = "La cresta brillante a V=0 (ZBCP ≈ 2e²/h) marca la fase topológica")
	vline!([-2, 2]; c = :white, ls = :dash, lw = 0.8, label = "")
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[compat]
Plots = "~1.41.6"
PlutoUI = "~0.7.83"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.2"
manifest_format = "2.0"
project_hash = "01875d4e5dc03778d4130a61534b5953cc3717b7"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "6c3913f4e9bdf6ba3c08041a446fb1332716cbc2"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.0"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitFlags]]
git-tree-sha1 = "bbe1079eecf9c9fbb52765193ad2bae27ae09bc8"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.10"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "1fa950ebc3e37eccd51c6a8fe1f92f7d86263522"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.7+0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "21d088c496ea22914fe80906eb5bce65755e5ec8"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.1"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e86f4a2805f7f19bec5129bc9150c38208e5dc23"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.4"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c307cd83373868391f3ac30b41530bc5d5d05d08"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.8.1+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "95ecf07c2eea562b5adbd0696af6db62c0f52560"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.5"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libva_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "cac41ca6b2d399adfc95e51240566f8a60a80806"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "8.1.0+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "70329abc09b886fd2c5d94ad2d9527639c421e3e"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.14.3+1"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "9e0fb9e54594c47f278d75063980e43066e26e20"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.1+1"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "44716a1a667cb867ee0e9ec8edc31c3e4aa5afdc"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.24"

    [deps.GR.extensions]
    IJuliaExt = "IJulia"

    [deps.GR.weakdeps]
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "be8a1b8065959e24fdc1b51402f39f3b6f0f6653"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.24+0"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "24f6def62397474a297bfcec22384101609142ed"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.3+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "51059d23c8bb67911a2e6fd5130229113735fc7e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.11.0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "f76f7560267b840e492180f9899b472f30b88450"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.6.0"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c0c9b76f3520863909825cbecdef58cd63de705a"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.5+0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "17b94ecafcfa45e8360a4fc9ca6b583b049e4e37"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.1.0+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cc3ad4faf30015a3e8094c9b5b7f19e85bdf2386"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.42.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "f04133fe05eff1667d2054c53d59f9122383fe05"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.2+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d620582b1f0cbe2c72dd1d5bd195a9ce73370ab1"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.42.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "77fe7779378a2331be7e86c64daaa2970bc2c1af"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "1.0.0"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "8785729fa736197687541f7053f6d8ab7fc44f92"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.10"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ff69a2b1330bcb730b9ac1ab7dd680176f5896b8"
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.1010+0"

[[deps.Measures]]
git-tree-sha1 = "b513cedd20d9c914783d8ad83d08120702bf2c77"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.5.20"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "NetworkOptions", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "1d1aaa7d449b58415f97d2839c318b70ffb525a0"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.6.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e2bb57a313a74b8104064b7efd01406c0a50d2ff"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.6.1+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.44.0+1"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58e5ed5e386e156bd93e86b305ebd21ac63d2d04"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.57.1+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "5d5e0a78e971354b1c7bff0655d11fdc1b0e12c8"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.4"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "e4a6721aa89e62e5d4217c0b21bd714263779dda"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.46.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "26ca162858917496748aad52bb5d3be4d26a228a"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.4"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "cb20a4eacda080e517e4deb9cfb6c7c518131265"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.41.6"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e189d0623e7ce9c37389bac17e80aac3b0302e75"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.83"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "edbeefc7a4889f528644251bdb5fc9ab5348bc2c"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "144895f6166994730ee7ff8113b981fc360638f1"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.10.2+2"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll", "Qt6Svg_jll"]
git-tree-sha1 = "d5b7dd0e226774cbd87e2790e34def09245c7eab"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.10.2+1"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "4d85eedf69d875982c46643f6b4f66919d7e157b"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.10.2+1"

[[deps.Qt6Svg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "81587ff5ff25a4e1115ce191e36285ede0334c9d"
uuid = "6de9746b-f93d-5813-b365-ba18ad4a9cf3"
version = "6.10.2+0"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "672c938b4b4e3e0169a07a5f227029d4905456f2"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.10.2+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "4f96c596b8c8258cc7d3b19797854d368f243ddc"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "c6f18e5a52a176a383f6f6c635e0f81feed1d6d4"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.11"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "82bee338d650aa515f31866c460cb7e3bcef90b8"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.8.2"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsStaticArraysCoreExt = ["StaticArraysCore"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "96478df35bbc2f3e1e791bc7a3d0eeee559e60e9"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.24.0+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b29c22e245d092b8b4e8d3c09ad7baa586d9f573"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.3+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3ea76ee3f4facd7a64684f9af25310825ee3668"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.2+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "9c7ad99c629a44f81e7799eb05ec2746abb5d588"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.6+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "808090ede1d41644447dd5cbafced4731c56bd2f"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.13+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c74ca84bbabc18c4547014765d194ff0b4dc9da"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.4+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "a376af5c7ae60d29825164db40787f15c80c7c54"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.3+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "0ba01bc7396896a4ace8aab67db31403c71628f4"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.7+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c174ef70c96c76f4c3f4d3cfbe09d018bcd1b53"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.6+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libpciaccess_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "58972370b81423fc546c56a60ed1a009450177c3"
uuid = "a65dc6b1-eb27-53a1-bb3e-dea574b5389e"
version = "0.19.0+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "ed756a03e95fff88d8f738ebc2849431bdd4fd1a"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.2.0+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "9750dc53819eba4e9a20be42349a6d3b86c7cdf8"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.6+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f4fc02e384b74418679983a97385644b67e1263b"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll"]
git-tree-sha1 = "68da27247e7d8d8dafd1fcf0c3654ad6506f5f97"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "44ec54b0e2acd408b0fb361e1e9244c60c9c3dd4"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "5b0263b6d080716a02544c55fdff2c8d7f9a16a0"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.10+0"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f233c83cad1fa0e70b7771e0e21b061a116f2763"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.2+0"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "ed349d26affcacafbc7fc2941ace1fb98f71e715"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.47.0+1"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c3b0e6196d50eab0c5ed34021aaa0bb463489510"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.14+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "850b06095ee71f0135d644ffd8a52850699581ed"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.3+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libdrm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "63aac0bcb0b582e11bad965cef4a689905456c03"
uuid = "8e53e030-5e6c-5a89-a30b-be5b7263a166"
version = "2.4.125+1"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56d643b57b188d30cccc25e331d416d3d358e557"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.13.4+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "91d05d7f4a9f67205bd6cf395e488009fe85b499"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.28.1+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e51150d5ab85cee6fc36726850f0e627ad2e4aba"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.58+0"

[[deps.libva_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "libdrm_jll"]
git-tree-sha1 = "7dbf96baae3310fe2fa0df0ccbb3c6288d5816c9"
uuid = "9a156e7d-b971-5f62-b2c9-67348b8fb97c"
version = "2.23.0+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b4d631fd51f2e9cdd93724ae25b2efc198b059b1"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "a1fc6507a40bf504527d0d4067d718f8e179b2b8"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.13.0+0"
"""

# ╔═╡ Cell order:
# ╟─fef64cb4-5d9f-11f1-86d6-b9cf870cb680
# ╟─fef755fa-5d9f-11f1-87ad-230064c776e7
# ╠═fef75640-5d9f-11f1-9745-0788a5745753
# ╠═fef7564a-5d9f-11f1-881a-4f94b82b21c4
# ╟─fef7567c-5d9f-11f1-9129-7730975337cb
# ╠═fef75690-5d9f-11f1-af3c-391647f08a24
# ╠═fef75698-5d9f-11f1-8fae-a3c138b5a123
# ╟─fef756a4-5d9f-11f1-961f-cb221fc79870
# ╠═fef756ae-5d9f-11f1-bc09-f5d1ea91bbeb
# ╟─fef756b8-5d9f-11f1-a111-279fef9838e0
# ╠═fef756cc-5d9f-11f1-9a43-819a539846c1
# ╟─fef756d6-5d9f-11f1-aaaf-778e2741c1a3
# ╠═fef756e0-5d9f-11f1-8626-d5fcfc82545e
# ╠═fef756ea-5d9f-11f1-99ae-6921d3e2e4f9
# ╟─fef7570a-5d9f-11f1-b08c-0b41129fc141
# ╠═fef75712-5d9f-11f1-880f-615c8eb0e4fc
# ╠═fef75726-5d9f-11f1-9dbc-0be9347a1a2b
# ╟─fef75730-5d9f-11f1-bd62-3de0a8163c23
# ╠═fef7573c-5d9f-11f1-ad24-93beaebf9fb0
# ╟─fef75744-5d9f-11f1-b376-f32645221222
# ╠═fef7574e-5d9f-11f1-a49e-63a88d4af75e
# ╠═fef75762-5d9f-11f1-a639-1b0db8fa8e2d
# ╟─fef75776-5d9f-11f1-91f8-ed0c37744c5d
# ╠═fef75780-5d9f-11f1-927f-e3cc3b6416a1
# ╟─fef75794-5d9f-11f1-8289-89f9d7298daa
# ╠═fef7579c-5d9f-11f1-9e8e-2ffee7f512e5
# ╠═fef757a8-5d9f-11f1-a24c-fdc49fe87dfa
# ╟─fef757bc-5d9f-11f1-b730-3b4c20c4c8e9
# ╠═fef757ce-5d9f-11f1-89d3-b5623a62ff4e
# ╠═fef757da-5d9f-11f1-aa5c-279ebcdd70e7
# ╠═fef757e4-5d9f-11f1-851d-a54019eca6b4
# ╠═fef757ee-5d9f-11f1-81bf-7501dcb720ad
# ╟─45be7096-5da4-11f1-9944-1913e66db732
# ╠═45bf77fc-5da4-11f1-9426-33509b8395ba
# ╠═45bf784c-5da4-11f1-8a05-7de32ed803cb
# ╟─45bf786a-5da4-11f1-8ca0-09dc71bb9675
# ╠═45bf7892-5da4-11f1-ab20-1f27e393a807
# ╟─fef757f8-5d9f-11f1-b13f-dbc40b051232
# ╠═fef75800-5d9f-11f1-8e03-4f7a65ada5ea
# ╠═fef75816-5d9f-11f1-95b6-df45280bd3b9
# ╟─fef75820-5d9f-11f1-95c3-91c9618303a1
# ╠═fef7582a-5d9f-11f1-83a7-41107188065b
# ╟─fef75834-5d9f-11f1-a38d-41975b7dcd99
# ╟─fef7583e-5d9f-11f1-b70b-6199465e6f42
# ╠═fef75852-5d9f-11f1-a9af-b9ea690da007
# ╠═fef7585c-5d9f-11f1-a52d-833c5e4f6150
# ╠═fef75866-5d9f-11f1-88fe-d77f4212fec6
# ╠═fef75872-5d9f-11f1-824f-13213e62ca60
# ╟─45d08718-5da4-11f1-8648-b95a4c5a96da
# ╠═45d087ae-5da4-11f1-b5d4-d75b6fb9092d
# ╠═45d087d8-5da4-11f1-92de-6926d279b7df
# ╟─fef7587a-5d9f-11f1-927f-b1b22b46a632
# ╠═fef7588e-5d9f-11f1-893a-1f654dc044e2
# ╠═fef75898-5d9f-11f1-81bb-cd27d653ddff
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
