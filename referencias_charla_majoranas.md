# Referencias — Qubits topológicos y los chips Majorana de Microsoft

**Charla UNIR · Máster en Computación Cuántica**
*Bibliografía de apoyo para estudiantes. Referencias verificadas a junio de 2026.*

> **Cómo leer esta lista.** Cada referencia indica revista, volumen/página y año, y el
> identificador **arXiv** (versión en acceso abierto) cuando existe. Las entradas marcadas
> con **★ REVIEW** son artículos de revisión o pedagógicos: el mejor punto de partida si
> quieres entender un concepto en profundidad. Para los preprints recientes de Microsoft se
> indica explícitamente si están **revisados por pares** o son **preprint**.

---

## ★ Para empezar (las mejores lecturas introductorias)

1. **M. A. Nielsen & I. L. Chuang**, *Quantum Computation and Quantum Information*, Cambridge University Press (2010). — El libro de texto estándar (qubits, puertas, corrección de errores).
2. **J. Preskill**, *"Quantum Computing in the NISQ era and beyond"*, **Quantum 2, 79 (2018)** — arXiv:1801.00862. ★ REVIEW del estado de la computación cuántica en su momento.
3. **R. Aguado**, *"Majorana quasiparticles in condensed matter"*, **Riv. Nuovo Cimento 40, 523 (2017)** — arXiv:1711.00011. ★ REVIEW pedagógica y accesible sobre Majoranas (muy recomendable).
4. **E. Prada *et al.***, *"From Andreev to Majorana bound states in hybrid superconductor–semiconductor nanowires"*, **Nat. Rev. Phys. 2, 575 (2020)** — arXiv:1911.04512. ★ REVIEW clave sobre por qué distinguir un Majorana real de un impostor trivial es tan difícil.

---

# BLOQUE 1 · Fundamentos: qubits lógicos, corrección de errores y hardware

### Conceptos y reviews
- **M. A. Nielsen & I. L. Chuang**, *Quantum Computation and Quantum Information*, Cambridge Univ. Press (2010). — Texto de referencia.
- **J. Preskill**, *"Quantum Computing in the NISQ era and beyond"*, **Quantum 2, 79 (2018)** — arXiv:1801.00862. ★ REVIEW.
- **B. M. Terhal**, *"Quantum error correction for quantum memories"*, **Rev. Mod. Phys. 87, 307 (2015)** — arXiv:1302.3428. ★ REVIEW de corrección de errores.
- **A. G. Fowler, M. Mariantoni, J. M. Martinis, A. N. Cleland**, *"Surface codes: Towards practical large-scale quantum computation"*, **Phys. Rev. A 86, 032324 (2012)** — arXiv:1208.0928. ★ REVIEW canónica del *surface code* y del coste (overhead).

### Origen de la corrección de errores cuántica
- **P. W. Shor**, *"Scheme for reducing decoherence in quantum computer memory"*, **Phys. Rev. A 52, R2493 (1995)**. — Primer código corrector cuántico.
- **A. M. Steane**, *"Error Correcting Codes in Quantum Theory"*, **Phys. Rev. Lett. 77, 793 (1996)**.
- **E. Dennis, A. Kitaev, A. Landahl, J. Preskill**, *"Topological quantum memory"*, **J. Math. Phys. 43, 4452 (2002)** — arXiv:quant-ph/0110143. — Origen del *surface code*.
- **W. K. Wootters & W. H. Zurek**, *"A single quantum cannot be cloned"*, **Nature 299, 802 (1982)**. — Teorema de no-clonación.

### Supremacía / ventaja cuántica
- **J. Preskill**, *"Quantum computing and the entanglement frontier"*, arXiv:1203.5813 (2012). — Acuña "quantum supremacy".
- **F. Arute *et al.* (Google)**, *"Quantum supremacy using a programmable superconducting processor"*, **Nature 574, 505 (2019)**. — Sycamore (RCS).
- **Google Quantum AI**, *"Observation of constructive interference at the edge of quantum ergodicity"* (algoritmo *Quantum Echoes*, OTOC), **Nature (2025)**, DOI 10.1038/s41586-025-09526-6. — Primera **ventaja cuántica verificable** (≈13.000×), octubre de 2025.

### Hitos de corrección de errores (experimento)
- **Google Quantum AI**, *"Quantum error correction below the surface code threshold"*, **Nature 638, 920–926 (2025)** — arXiv:2408.13687. — Willow: primer *below-threshold* (Λ = 2,14), diciembre de 2024.
- **S. Bravyi *et al.* (IBM)**, *"High-threshold and low-overhead fault-tolerant quantum memory"*, **Nature 627, 778 (2024)** — arXiv:2308.07915. — Códigos qLDPC (overhead reducido).

### El "blanco que se mueve": coste de romper RSA
- **C. Gidney & M. Ekerå**, *"How to factor 2048 bit RSA integers in 8 hours using 20 million noisy qubits"*, **Quantum 5, 433 (2021)** — arXiv:1905.09749.
- **C. Gidney**, *"How to factor 2048 bit RSA integers with less than a million noisy qubits"*, arXiv:2505.15917 (2025). — La estimación cae de 20M a <1M en 6 años.

---

# BLOQUE 2 · Modos de Majorana: teoría, materiales y firmas

### Origen y modelo de Kitaev
- **E. Majorana**, *"Teoria simmetrica dell'elettrone e del positrone"*, **Nuovo Cimento 14, 171 (1937)**. — La partícula que es su propia antipartícula.
- **A. Yu. Kitaev**, *"Unpaired Majorana fermions in quantum wires"*, **Phys.-Usp. 44, 131 (2001)** — arXiv:cond-mat/0010440. — El modelo de juguete (cadena 1D).

### El nanohilo semiconductor como realización del modelo de Kitaev
*(las dos referencias canónicas, siempre citadas juntas)*
- **R. M. Lutchyn, J. D. Sau, S. Das Sarma**, *"Majorana Fermions and a Topological Phase Transition in Semiconductor-Superconductor Heterostructures"*, **Phys. Rev. Lett. 105, 077001 (2010)** — arXiv:1002.4033.
- **Y. Oreg, G. Refael, F. von Oppen**, *"Helical Liquids and Majorana Bound States in Quantum Wires"*, **Phys. Rev. Lett. 105, 177002 (2010)** — arXiv:1003.1145.

*Precursores:*
- **J. D. Sau, R. M. Lutchyn, S. Tewari, S. Das Sarma**, *"Generic New Platform for Topological Quantum Computation Using Semiconductor Heterostructures"*, **Phys. Rev. Lett. 104, 040502 (2010)** — arXiv:0907.2239.
- **L. Fu & C. L. Kane**, *"Superconducting Proximity Effect and Majorana Fermions at the Surface of a Topological Insulator"*, **Phys. Rev. Lett. 100, 096407 (2008)** — arXiv:0707.1692.

### Reviews de Majoranas y superconductividad topológica  ★
- **R. Aguado**, *"Majorana quasiparticles in condensed matter"*, **Riv. Nuovo Cimento 40, 523 (2017)** — arXiv:1711.00011. ★ La más pedagógica; ideal para empezar.
- **J. Alicea**, *"New directions in the pursuit of Majorana fermions in solid state systems"*, **Rep. Prog. Phys. 75, 076501 (2012)** — arXiv:1202.1293. ★ REVIEW clásica.
- **M. Leijnse & K. Flensberg**, *"Introduction to topological superconductivity and Majorana fermions"*, **Semicond. Sci. Technol. 27, 124003 (2012)** — arXiv:1206.1736. ★ Introducción breve y clara.
- **C. W. J. Beenakker**, *"Search for Majorana fermions in superconductors"*, **Annu. Rev. Condens. Matter Phys. 4, 113 (2013)** — arXiv:1112.1950. ★ REVIEW.
- **M. Sato & Y. Ando**, *"Topological superconductors: a review"*, **Rep. Prog. Phys. 80, 076501 (2017)** — arXiv:1608.03395. ★ REVIEW más amplia.
- **R. M. Lutchyn *et al.***, *"Majorana zero modes in superconductor–semiconductor heterostructures"*, **Nat. Rev. Mater. 3, 52 (2018)** — arXiv:1707.04899. ★ REVIEW específica de la plataforma de nanohilos.

### Computación topológica y anyones no-Abelianos  ★
- **C. Nayak, S. H. Simon, A. Stern, M. Freedman, S. Das Sarma**, *"Non-Abelian anyons and topological quantum computation"*, **Rev. Mod. Phys. 80, 1083 (2008)** — arXiv:0707.1889. ★ REVIEW de referencia (braiding, fusión, universalidad).

### La firma: pico a voltaje cero (ZBCP) y reflexión de Andreev
- **K. T. Law, P. A. Lee, K. T. Ng**, *"Majorana Fermion Induced Resonant Andreev Reflection"*, **Phys. Rev. Lett. 103, 237001 (2009)** — arXiv:0907.1909. — Predice el pico cuantizado a 2e²/h.

### Experimentos pioneros (firmas en nanohilos)
- **V. Mourik *et al.***, *"Signatures of Majorana fermions in hybrid superconductor-semiconductor nanowire devices"*, **Science 336, 1003 (2012)**. — Primer experimento (nanohilo de InSb).
- **A. Das *et al.***, *"Zero-bias peaks and splitting in an Al–InAs nanowire topological superconductor…"*, **Nat. Phys. 8, 887 (2012)**. — Plataforma **InAs/Al**.
- **M. T. Deng *et al.***, *"Majorana bound state in a coupled quantum-dot hybrid-nanowire system"*, **Science 354, 1557 (2016)** — arXiv:1612.07989.

### El problema de los falsos positivos (estados triviales / quasi-Majorana)
*(esencial para entender por qué la cuantización del pico no basta)*
- **G. Kells, D. Meidan, P. W. Brouwer**, *"Near-zero-energy end states in topologically trivial spin-orbit coupled superconducting nanowires…"*, **Phys. Rev. B 86, 100503(R) (2012)**.
- **C. Vuik, B. Nijholt, A. R. Akhmerov, M. Wimmer**, *"Reproducing topological properties with quasi-Majorana states"*, **SciPost Phys. 7, 061 (2019)** — arXiv:1806.02801.
- **H. Pan & S. Das Sarma**, *"Physical mechanisms for zero-bias conductance peaks in Majorana nanowires"*, **Phys. Rev. Research 2, 013377 (2020)** — arXiv:1912.07747.
- **E. Prada *et al.***, Nat. Rev. Phys. 2, 575 (2020) — arXiv:1911.04512. ★ REVIEW (citada arriba): Andreev vs Majorana.

---

# BLOQUE 3 · Microsoft, el tetrón y los chips Majorana 1/2

### Computación por medidas (sin braiding) y la arquitectura del tetrón
- **P. Bonderson, M. Freedman, C. Nayak**, *"Measurement-Only Topological Quantum Computation"*, **Phys. Rev. Lett. 101, 010501 (2008)** — arXiv:0802.0279. — Base del esquema "operar midiendo".
- **T. Karzig *et al.***, *"Scalable designs for quasiparticle-poisoning-protected topological quantum computation with Majorana zero modes"*, **Phys. Rev. B 95, 235305 (2017)** — arXiv:1610.05289. — **Propuesta original del tetrón** (acuña los términos "tetrón" y "hexón").
- **S. Plugge, A. Rasmussen, R. Egger, K. Flensberg**, *"Majorana box qubits"*, **New J. Phys. 19, 012001 (2017)** — arXiv:1609.01697. — Propuesta equivalente y contemporánea.

### El Topological Gap Protocol (TGP) y el debate
- **D. I. Pikulin *et al.***, *"Protocol to identify a topological superconducting phase in a three-terminal device"*, arXiv:2103.12217 (2021). — Propuesta del TGP.
- **Microsoft Quantum (M. Aghaee *et al.*)**, *"InAs-Al hybrid devices passing the topological gap protocol"*, **Phys. Rev. B 107, 245423 (2023)** — arXiv:2207.02472. *(Revisado por pares; publicado el 21 de junio de 2023.)*
- **H. F. Legg**, *"Comment on 'InAs-Al hybrid devices passing the topological gap protocol'…"*, arXiv:2502.19560 (2025). — Crítica: riesgo de falsos positivos según la elección de parámetros.

### Los chips Majorana 1 y 2 (papers primarios)
- **Microsoft Azure Quantum (M. Aghaee *et al.*)**, *"Interferometric single-shot parity measurement in InAs–Al hybrid devices"*, **Nature 638, 651–655 (2025)** — DOI 10.1038/s41586-024-08445-2. *(★ **Revisado por pares**. Publicado el **19 de febrero de 2025**, junto al anuncio de Majorana 1. Lleva una **nota editorial** advirtiendo de que los resultados no constituyen evidencia de modos de Majorana. Precursor: arXiv:2401.09549.)*
- **Microsoft Quantum**, *"Distinct Lifetimes for X and Z Loop Measurements in a Majorana Tetron Device"*, arXiv:2507.08795 (**julio de 2025**). *(Preprint. Tetrón con los bucles X y Z; τ_X ≈ 14,5 µs, τ_Z ≈ 12,4 ms.)*
- **Microsoft Quantum**, *"20 Second Parity Lifetime in an InAs–Pb Tetron Device"*, arXiv:2606.03884 (**2 de junio de 2026**). *(Preprint. Majorana 2: plomo en vez de aluminio, vida de paridad Z ~20 s.)*

### Contexto histórico: la retractación de 2018
- **H. Zhang *et al.***, *"Quantized Majorana conductance"*, **Nature 556, 74 (2018)** — **RETRACTADO en 2021** (*Nature 591, E30*).

### Hoja de ruta
- **Microsoft Quantum**, *"Roadmap to fault-tolerant quantum computation using topological qubit arrays"*, arXiv:2502.12252 (febrero de 2025). — Las cuatro generaciones de dispositivo hacia el prototipo tolerante a fallos.

---

## Nota sobre el estado de revisión (importante para evaluar la evidencia)

| Trabajo | Estado |
|---|---|
| Nature 638, 651 (2025) — Majorana 1 (lectura single-shot) | 🟢 **Revisado por pares** (con nota editorial) |
| PRB 107, 245423 (2023) — TGP | 🟢 **Revisado por pares** |
| arXiv:2507.08795 (jul 2025) — tetrón X/Z | 🟠 **Preprint** |
| arXiv:2606.03884 (jun 2026) — Majorana 2 | 🟠 **Preprint** |
| arXiv:2502.19560 (2025) — crítica de Legg | 🟠 **Preprint (comment)** |

> A junio de 2026 **no hay consenso** en la comunidad sobre si estos dispositivos albergan
> genuinamente modos de Majorana. La distinción entre **promesa estructural** (real y enorme
> si funciona) y **evidencia demostrada** (en disputa) es un aspecto crucial de la computación cuántica topológica basada en Majoranas.
