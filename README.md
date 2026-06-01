# Modos de Majorana en la cadena de Kitaev

Notebook interactivo (en **Julia**, con [Pluto.jl](https://plutojl.org)) que construye paso
a paso la física de los modos de Majorana en el modelo de Kitaev 1D: espectro y modo de
energía cero, base de Majorana, visualización de los modos de borde, función espectral /
LDOS, diagramas de fase, invariantes topológicos del bulk (Pfaffiano `ℤ₂`, número de
winding y fase de Zak `ℤ`) y el pico de conductancia a voltaje cero (ZBCP), incluyendo un
cálculo de transporte cuantizado (`2e²/h`) con la reflexión de Andreev **implementado en
Julia puro**, sin dependencias externas.

- **Notebook:** [`majoranas_in_1D.jl`](majoranas_in_1D.jl) — notebook reactivo de Pluto.

El notebook es **reactivo**: mueve los *sliders* y todas las figuras se recalculan al
instante.

---

## Descarga y última versión

La **versión más reciente** de este notebook está siempre en GitHub:

👉 **https://github.com/mas9/intro-modos-majorana**

> Si te han pasado esta carpeta en un ZIP, puede quedar desactualizada. En el repositorio
> encontrarás siempre la última versión.

Tienes dos formas de conseguirlo:

- **Descargar el ZIP** (no necesitas git): abre el repositorio, pulsa el botón verde
  **`Code` → `Download ZIP`** y descomprime la carpeta.
- **Clonar con git** (recomendado para actualizarlo cómodamente):
  ```bash
  git clone https://github.com/mas9/intro-modos-majorana.git
  cd intro-modos-majorana
  ```
  Para traerte los cambios publicados más adelante, desde la carpeta del proyecto:
  ```bash
  git pull
  ```

---

## Requisitos

- **Julia ≥ 1.10** (probado con 1.12). Recomendado instalarla con
  [`juliaup`](https://github.com/JuliaLang/juliaup):
  ```bash
  curl -fsSL https://install.julialang.org | sh     # macOS / Linux
  # Windows:  winget install julia -s msstore
  ```

No hace falta instalar nada más a mano: **Pluto descarga e instala automáticamente todos
los paquetes en las versiones exactas** que necesita el notebook. Esas versiones están
**fijadas dentro del propio `.jl`** (Pluto guarda un `Project.toml` y un `Manifest.toml`
embebidos al final del archivo), así que el notebook funciona igual en cualquier ordenador.

> En el **primer arranque** Pluto descargará y precompilará los paquetes (Plots, PlutoUI…).
> Puede tardar varios minutos; las siguientes veces es instantáneo.

---

## Cómo ejecutarlo

### Opción recomendada (la más sencilla)

1. Instala Pluto una sola vez (en el entorno global de Julia):
   ```bash
   julia -e 'using Pkg; Pkg.add("Pluto")'
   ```
2. Arranca Pluto:
   ```bash
   julia -e 'using Pluto; Pluto.run()'
   ```
   Se abrirá el navegador. En la pantalla de inicio, en *"Open a notebook"*, pega la ruta a
   `majoranas_in_1D.jl` (o navega hasta él) y ábrelo. Pluto instalará automáticamente las
   dependencias embebidas y ejecutará el notebook.

### Opción con el entorno de este proyecto

Este repositorio incluye un `Project.toml` y un `Manifest.toml` que fijan una versión
concreta de Pluto y los paquetes para *lanzar* el notebook. Desde la carpeta del proyecto:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'      # instala el entorno (una vez)
julia --project=. -e 'using Pluto; Pluto.run()'          # arranca Pluto
```

y abre `majoranas_in_1D.jl` desde la interfaz de Pluto.

---

## Estructura del repositorio

| Archivo | Descripción |
|---|---|
| `majoranas_in_1D.jl`   | **Notebook de Pluto** (en Julia). Lleva su entorno embebido. |
| `Project.toml`         | Dependencias para lanzar Pluto (`Pluto`, `Plots`, `PlutoUI`) con límites de versión. |
| `Manifest.toml`        | Versiones exactas resueltas del entorno anterior (reproducibilidad del *lanzador*). |

---

## Contenido del notebook

0. Preparativos
1. El Hamiltoniano de Bogoliubov-de Gennes (BdG)
2. Espectro: la aparición del modo de energía cero
3. De fermiones a Majoranas (base de Majorana)
4. Visualización de los modos de Majorana en los extremos
5. Función espectral y densidad local de estados (LDOS)
6. Regímenes de `t`, `μ` y `Δ`: diagrama de fases
7. Invariantes del bulk: el Pfaffiano (`ℤ₂`), el número de winding y la fase de Zak (`ℤ`)
8. El pico de conductancia a voltaje cero (ZBCP) — proxy espectral
9. Resumen
10. Conductancia cuantizada `2e²/h` y reflexión de Andreev (transporte en Julia puro)

---

## Solución de problemas

- **"No me funcionan los *sliders*":** asegúrate de abrir el `.jl` **dentro de Pluto**, no
  como script. Los controles `@bind` solo son interactivos en Pluto.
- **Primer arranque muy lento:** es normal (descarga + precompilación). Espera a que
  termine; no cierres la pestaña.
- **Quiero reiniciar el entorno embebido:** en Pluto, menú de paquetes (icono de cajas
  arriba a la derecha) → *"Update packages"* o borra los bloques
  `PLUTO_PROJECT_TOML_CONTENTS` / `PLUTO_MANIFEST_TOML_CONTENTS` del final del `.jl` y
  vuelve a abrirlo.
