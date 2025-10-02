Overview

This repository provides a set of Python scripts for simulating the combustion behavior of hydrogen-enriched fuels (methane, methanol, propane) using Cantera, combined with thermodynamic validation. The simulations assess:

Heat release rates and total heat released

Adiabatic flame temperatures

NOx and CO₂ emissions

Ignition behavior across residence times

Influence of equivalence ratio (ϕ) and hydrogen blending levels

The scripts illustrate the transition from conventional hydrocarbon fuels to hydrogen-enriched blends, supporting practical strategies for decarbonization of engines, turbines, and CHP systems.

Features

Hydrogen blending analysis

Methanol–Hydrogen, Methane–Hydrogen, Propane–Hydrogen mixtures

Hydrogen fractions from 10% to 50% (by mole)

Equivalence ratios ϕ = 0.75, 1.0, 1.5

Performance metrics

Total heat release (MJ/m³ and kJ)

Peak NOx and CO₂ emissions (ppm)

Final combustion temperature (K)

Ignition timing across short residence times

Validation

Comparison between thermodynamic HHV-based models and Cantera equilibrium/kinetics

Heat release trends match within 10–15%

Visualization

Heat maps (H₂ % vs. ϕ) for heat release, NOx, CO₂

Temperature rise vs. residence time (log scale)

Heat release rate vs. time

LHV and heat per kg air for blended fuels

Requirements

Python 3.8+

Cantera
 ≥ 3.0

NumPy, Matplotlib, Pandas

Install dependencies:

pip install cantera numpy matplotlib pandas

Usage

Each script focuses on a different fuel blend and output:

methanol_h2.py

Simulates CH₃OH–H₂ blends

Outputs heat release, NOx emissions, and temperature histories

methane_h2.py

Simulates CH₄–H₂ blends

Adds CO₂ emission tracking

propane_h2.py

Simulates C₃H₈–H₂ blends

Provides heat release and NOx heatmaps

thermo_validation.py

Compares Cantera equilibrium with thermodynamic HHV-based estimates

Calculates mixture LHV and energy release per kg of air

ignition_dynamics.py

Evaluates ignition behavior for CH₄–H₂ blends

Varies spark temperature (1100–1200 K) and residence time (1–5 ms)

Run any script with:

python methane_h2.py

Outputs

Printed results: heat (kJ), final T (K), NOx and CO₂ (ppm) per blend ratio

Figures:

Heat released vs H₂%

NOx/CO₂ vs H₂%

Final T vs H₂%

Heat maps (ϕ vs. H₂%)

Time-resolved curves for temperature and heat release

Applications

Provides validated insights for:

Hydrogen co-firing in gas turbines

Hydrogen blending in SI and CI engines

CHP integration with hydrogen engines

Serves as a fast-track computational roadmap for experimental calibration, retrofits, and emissions compliance.

Citation

If you use these codes or results, please cite:

Abbass, A. Hydrogen Blending and Engine/Turbine Transition Pathways for Decarbonized Power and CHP Systems. Mississippi State University, 2025.
