# Ten Traveltime Equations in Seismics: Complete Derivations and Applications

[![Python](https://img.shields.io/badge/Language-Python-3776AB?style=flat&logo=python&logoColor=white)](https://www.python.org/)
[![Jupyter Notebook](https://img.shields.io/badge/Tools-Jupyter%20Notebook-F37626?style=flat&logo=jupyter&logoColor=white)](https://jupyter.org/)
[![Journal](https://img.shields.io/badge/Journal-BrJG__Submit-blue.svg)](https://sbgf.org.br/revista/)

Traveltime calculation is of fundamental importance in seismic processing, interpretation, imaging, and tomography. This repository consolidates and expands a pedagogical review containing **ten essential traveltime equations** commonly introduced in undergraduate geophysics. 

Rather than starting from abstract eikonal solver formulations, this project uses an intentionally pedagogical approach. Every derivation emphasizes geometry, trigonometry, image-point constructions, and physical intuition, connecting simple analytical proofs to practical processing workflows.

---

## 📋 Project Structure & Covered Equations

The interactive Jupyter Notebook maps directly to the paper's mathematical derivations, exploring a two-layer, homogeneous, isotropic model governed by velocities ($v_1, v_2$), interface depth ($z$), and acquisition offset ($x$):

| Section | Seismic Configuration | Mathematical & Geometric Core | Practical Application |
| :--- | :--- | :--- | :--- |
| **1** | **Direct Waves** | Linear propagation: $T = \frac{x}{v_1}$ | First-arrival analysis |
| **2** | **Primary Reflections** | Pythagorean theorem & reflection hyperbola | Normal Moveout (NMO) |
| **3** | **First-Order Multiples** | Extended image-point construction | Multiple prediction & attenuation |
| **4** | **Dipping Plane Reflector** | Law of Cosines & reflector dip ($\theta$) | Asymmetric moveout analysis |
| **5** | **Dipping First-Order Multiples** | Multi-stage image-point mapping | Complex multiple attenuation |
| **6** | **Converted $P\text{-}SV$ / $SV\text{-}P$ Waves** | Segment-by-segment mode coupling | Multicomponent processing |
| **7** | **Refracted Waves (Horizontal)** | Critical angle ($\theta_c$) and head-wave path | Shallow refraction statics |
| **8** | **Refracted Waves (Dipping)** | Updip/downdip asymmetry relations | Refractor dip & velocity estimation |
| **9** | **Diffractions (Outside Intercept)** | Double-Square-Root (DSR) equation | Kirchhoff migration |
| **10**| **Diffractions (Inside Intercept)** | Midpoint-offset coordinate formulation | Diffraction imaging |
| **Bonus**| **1D Heterogeneous Medium** | Linear velocity gradient $v(z) = v_0 + az$ | Continuous ray tracing |

---

## 🛠️ Features of the Jupyter Notebook

The code matches the pedagogical spirit of the manuscript:
* **Interactive Sliders:** Modify parameters like layer velocities ($v_1, v_2$), interface dips ($\alpha, \theta$), or velocity gradients ($a$) using `ipywidgets` and observe real-time deformation of traveltime curves.
* **Hyperbolic Moveout Overlays:** Visually compare how dip or wave conversion deviates the responses from standard horizontal reflection hyperbolas.
* **Continuum Limit Verification:** Check numerically how discretization turns into analytical continuous solutions in the 1D linear gradient model.

---

## 🚀 Getting Started

### Prerequisites
To run the notebooks locally, you need Python 3.x and the following core scientific packages:
* `numpy`
* `matplotlib`
* `ipywidgets` (for the interactive visualization components)

### Installation & Execution
1. Clone this repository:
   ```bash
   git clone [https://github.com/annie-gabrielle/10_Seismic_Equations.git](https://github.com/annie-gabrielle/10_Seismic_Equations.git)
