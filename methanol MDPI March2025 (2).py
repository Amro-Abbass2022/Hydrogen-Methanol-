#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# Setup
gas = ct.Solution('gri30.yaml', transport_model=None)
equivalence_ratios = [0.75, 1.0, 1.5]
h2_fractions = np.linspace(0.1, 0.5, 5)

results = {phi: {'H2%': [], 'Heat_kJ': [], 'NOx_ppm': [], 'Final_Temp': [],
                 'TimeHistories': {}} for phi in equivalence_ratios}

# Run simulations
for phi in equivalence_ratios:
    print(f"\n--- φ = {phi} ---")
    for h2_frac in h2_fractions:
        ch3oh_frac = 1 - h2_frac
        fuel_mixture = f'CH3OH:{ch3oh_frac}, H2:{h2_frac}'
        gas.TP = 300.0, ct.one_atm
        gas.set_equivalence_ratio(phi, fuel=fuel_mixture, oxidizer='O2:1.0, N2:3.76')

        inlet = ct.Reservoir(gas)
        gas.equilibrate('HP')
        combustor = ct.IdealGasReactor(gas)
        combustor.volume = 1.0
        exhaust = ct.Reservoir(gas)

        def mdot(t): return combustor.mass / residence_time
        inlet_mfc = ct.MassFlowController(inlet, combustor, mdot=mdot)
        outlet_mfc = ct.PressureController(combustor, exhaust, primary=inlet_mfc, K=0.01)

        sim = ct.ReactorNet([combustor])
        residence_time = 0.1

        times, heat_rates, nox_vals, temps = [], [], [], []

        while combustor.T > 500:
            sim.initial_time = 0.0
            sim.advance_to_steady_state()
            times.append(residence_time)
            heat_rates.append(combustor.thermo.heat_release_rate)
            temps.append(combustor.T)
            NOx_ppm = (combustor.thermo['NO'].X[0] + combustor.thermo['NO2'].X[0]) * 1e6
            nox_vals.append(NOx_ppm)
            residence_time *= 0.9

        Q_kJ = np.trapz(heat_rates, times) * combustor.volume / 1000
        NOx_peak = max(nox_vals)
        T_final = combustor.T
        H2_percent = round(h2_frac * 100, 1)

        results[phi]['H2%'].append(H2_percent)
        results[phi]['Heat_kJ'].append(Q_kJ)
        results[phi]['NOx_ppm'].append(NOx_peak)
        results[phi]['Final_Temp'].append(T_final)
        results[phi]['TimeHistories'][H2_percent] = {
            'time': times, 'temp': temps, 'heat': heat_rates
        }

        print(f"H₂ {H2_percent:.0f}% → Heat: {Q_kJ:.2f} kJ, Temp: {T_final:.1f} K, NOx: {NOx_peak:.1f} ppm")

# === Original Summary Plots ===

# Total Heat Released vs H₂%
plt.figure(figsize=(8, 5))
for phi in equivalence_ratios:
    plt.plot(results[phi]['H2%'], results[phi]['Heat_kJ'], marker='o', label=f'ϕ = {phi}')
plt.xlabel('Hydrogen %')
plt.ylabel('Total Heat Released (kJ)')
plt.title('Total Heat Released vs Hydrogen % (All φ)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# NOx vs H₂%
plt.figure(figsize=(8, 5))
for phi in equivalence_ratios:
    plt.plot(results[phi]['H2%'], results[phi]['NOx_ppm'], marker='s', label=f'ϕ = {phi}')
plt.xlabel('Hydrogen %')
plt.ylabel('NOx Emissions (ppm)')
plt.title('NOx Emissions vs Hydrogen % (All φ)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Final Temperature vs H₂%
plt.figure(figsize=(8, 5))
for phi in equivalence_ratios:
    plt.plot(results[phi]['H2%'], results[phi]['Final_Temp'], marker='^', label=f'ϕ = {phi}')
plt.xlabel('Hydrogen %')
plt.ylabel('Final Temperature (K)')
plt.title('Final Temperature vs Hydrogen % (All φ)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Heat vs Final Temperature
plt.figure(figsize=(8, 5))
for phi in equivalence_ratios:
    plt.plot(results[phi]['Final_Temp'], results[phi]['Heat_kJ'], marker='d', label=f'ϕ = {phi}')
plt.xlabel('Final Temperature (K)')
plt.ylabel('Total Heat Released (kJ)')
plt.title('Heat Released vs Final Temperature (per φ)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# === NEW: Temperature Rise and Heat Release vs Time ===

for phi in equivalence_ratios:
    # Temperature Rise
    plt.figure(figsize=(8, 5))
    for H2_percent, data in results[phi]['TimeHistories'].items():
        plt.plot(data['time'], data['temp'], marker='o', label=f"H₂ {H2_percent}%")
    plt.xscale('log')
    plt.xlabel("Residence Time [s] (log scale)")
    plt.ylabel("Temperature [K]")
    plt.title(f"Temperature Rise vs Time (ϕ = {phi})")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Heat Release Rate
    plt.figure(figsize=(8, 5))
    for H2_percent, data in results[phi]['TimeHistories'].items():
        plt.plot(data['time'], data['heat'], marker='s', label=f"H₂ {H2_percent}%")
    plt.xscale('log')
    plt.xlabel("Residence Time [s] (log scale)")
    plt.ylabel("Heat Release Rate [W/m³]")
    plt.title(f"Heat Release Rate vs Time (ϕ = {phi})")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


# In[2]:


import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# Setup
gas = ct.Solution('gri30.yaml', transport_model=None)
equivalence_ratios = [0.75, 1.0, 1.5]
h2_fractions = np.linspace(0.1, 0.5, 5)

results = {phi: {
    'H2%': [], 'Heat_kJ': [], 'NOx_ppm': [], 'CO2_ppm': [], 'Final_Temp': [],
    'TimeHistories': {}
} for phi in equivalence_ratios}

# Run simulations
for phi in equivalence_ratios:
    print(f"\n--- φ = {phi} ---")
    for h2_frac in h2_fractions:
        ch4_frac = 1 - h2_frac
        fuel_mixture = f'CH4:{ch4_frac}, H2:{h2_frac}'
        gas.TP = 300.0, ct.one_atm
        gas.set_equivalence_ratio(phi, fuel=fuel_mixture, oxidizer='O2:1.0, N2:3.76')

        inlet = ct.Reservoir(gas)
        gas.equilibrate('HP')
        combustor = ct.IdealGasReactor(gas)
        combustor.volume = 1.0
        exhaust = ct.Reservoir(gas)

        def mdot(t): return combustor.mass / residence_time
        inlet_mfc = ct.MassFlowController(inlet, combustor, mdot=mdot)
        outlet_mfc = ct.PressureController(combustor, exhaust, primary=inlet_mfc, K=0.01)

        sim = ct.ReactorNet([combustor])
        residence_time = 0.1

        times, heat_rates, nox_vals, temps, co2_vals = [], [], [], [], []

        while combustor.T > 500:
            sim.initial_time = 0.0
            sim.advance_to_steady_state()
            times.append(residence_time)
            heat_rates.append(combustor.thermo.heat_release_rate)
            temps.append(combustor.T)
            NOx_ppm = (combustor.thermo['NO'].X[0] + combustor.thermo['NO2'].X[0]) * 1e6
            CO2_ppm = combustor.thermo['CO2'].X[0] * 1e6
            nox_vals.append(NOx_ppm)
            co2_vals.append(CO2_ppm)
            residence_time *= 0.9

        Q_kJ = np.trapz(heat_rates, times) * combustor.volume / 1000
        NOx_peak = max(nox_vals)
        CO2_peak = max(co2_vals)
        T_final = combustor.T
        H2_percent = round(h2_frac * 100, 1)

        results[phi]['H2%'].append(H2_percent)
        results[phi]['Heat_kJ'].append(Q_kJ)
        results[phi]['NOx_ppm'].append(NOx_peak)
        results[phi]['CO2_ppm'].append(CO2_peak)
        results[phi]['Final_Temp'].append(T_final)
        results[phi]['TimeHistories'][H2_percent] = {
            'time': times, 'temp': temps, 'heat': heat_rates
        }

        print(f"H₂ {H2_percent:.0f}% → Heat: {Q_kJ:.2f} kJ, Temp: {T_final:.1f} K, NOx: {NOx_peak:.1f} ppm, CO₂: {CO2_peak:.1f} ppm")

# === Summary Plots ===

# Heat Released vs H₂%
plt.figure(figsize=(8, 5))
for phi in equivalence_ratios:
    plt.plot(results[phi]['H2%'], results[phi]['Heat_kJ'], marker='o', label=f'ϕ = {phi}')
plt.xlabel('Hydrogen %')
plt.ylabel('Total Heat Released (kJ)')
plt.title('Total Heat Released vs Hydrogen % (CH₄ + H₂)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# NOx vs H₂%
plt.figure(figsize=(8, 5))
for phi in equivalence_ratios:
    plt.plot(results[phi]['H2%'], results[phi]['NOx_ppm'], marker='s', label=f'ϕ = {phi}')
plt.xlabel('Hydrogen %')
plt.ylabel('NOx Emissions (ppm)')
plt.title('NOx Emissions vs Hydrogen % (CH₄ + H₂)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# CO₂ vs H₂%
plt.figure(figsize=(8, 5))
for phi in equivalence_ratios:
    plt.plot(results[phi]['H2%'], results[phi]['CO2_ppm'], marker='^', label=f'ϕ = {phi}')
plt.xlabel('Hydrogen %')
plt.ylabel('CO₂ Emissions (ppm)')
plt.title('CO₂ Emissions vs Hydrogen % (CH₄ + H₂)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Final Temperature vs H₂%
plt.figure(figsize=(8, 5))
for phi in equivalence_ratios:
    plt.plot(results[phi]['H2%'], results[phi]['Final_Temp'], marker='^', label=f'ϕ = {phi}')
plt.xlabel('Hydrogen %')
plt.ylabel('Final Temperature (K)')
plt.title('Final Temperature vs Hydrogen % (CH₄ + H₂)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# === Time-Resolved Plots: Temperature and Heat ===

for phi in equivalence_ratios:
    # Temperature Rise
    plt.figure(figsize=(8, 5))
    for H2_percent, data in results[phi]['TimeHistories'].items():
        plt.plot(data['time'], data['temp'], marker='o', label=f"H₂ {H2_percent}%")
    plt.xscale('log')
    plt.xlabel("Residence Time [s]")
    plt.ylabel("Temperature [K]")
    plt.title(f"Temperature vs Time (ϕ = {phi})")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Heat Release Rate
    plt.figure(figsize=(8, 5))
    for H2_percent, data in results[phi]['TimeHistories'].items():
        plt.plot(data['time'], data['heat'], marker='s', label=f"H₂ {H2_percent}%")
    plt.xscale('log')
    plt.xlabel("Residence Time [s]")
    plt.ylabel("Heat Release Rate [W/m³]")
    plt.title(f"Heat Release vs Time (ϕ = {phi})")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


# In[3]:


import numpy as np
import matplotlib.pyplot as plt

# Data from uploaded plots

hydrogen_percent = [10, 20, 30, 40, 50]
phis = [0.75, 1.0, 1.5]

# Methanol-Hydrogen Mixture (Heat Released, NOx)
methanol_heat = np.array([
    [-2470, -2490, -2535, -2640, -2720],
    [-2820, -2860, -2900, -2950, -3030],
    [-2530, -2570, -2610, -2700, -2750]
])
methanol_nox = np.array([
    [110, 115, 120, 135, 150],
    [860, 900, 950, 1010, 1100],
    [80, 82, 85, 87, 90]
])

# Methane-Hydrogen Mixture (Heat Released, NOx, CO2)
methane_heat = np.array([
    [-2430, -2470, -2540, -2600, -2670],
    [-2760, -2800, -2840, -2910, -2970],
    [-1710, -1810, -1910, -2020, -2140]
])
methane_nox = np.array([
    [75, 80, 85, 95, 105],
    [870, 910, 950, 1000, 1080],
    [45, 50, 60, 70, 85]
])
methane_co2 = np.array([
    [70200, 67800, 64800, 61200, 56800],
    [81200, 78000, 74400, 69800, 64600],
    [40200, 38400, 36200, 33600, 30600]
])

def plot_heatmap(data, title, ylabel, cmap='viridis'):
    plt.figure(figsize=(8, 5))
    X, Y = np.meshgrid(hydrogen_percent, phis)
    c = plt.contourf(X, Y, data, levels=20, cmap=cmap)
    plt.xlabel('Hydrogen %')
    plt.ylabel('Equivalence Ratio (ϕ)')
    plt.title(title)
    plt.colorbar(c, label=ylabel)
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.tight_layout()
    plt.show()

# === Heatmaps ===
plot_heatmap(methanol_heat, 'Heat Map for Methanol-Hydrogen (Heat Released)', 'Heat Released (kJ)')
plot_heatmap(methanol_nox, 'Heat Map for Methanol-Hydrogen (NOx Emissions)', 'NOx (ppm)', cmap='plasma')

plot_heatmap(methane_heat, 'Heat Map for Methane-Hydrogen (Heat Released)', 'Heat Released (kJ)')
plot_heatmap(methane_nox, 'Heat Map for Methane-Hydrogen (NOx Emissions)', 'NOx (ppm)', cmap='plasma')
plot_heatmap(methane_co2, 'Heat Map for Methane-Hydrogen (CO₂ Emissions)', 'CO₂ (ppm)', cmap='inferno')


# In[ ]:






# In[1]:


import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# Setup
gas = ct.Solution('gri30.yaml', transport_model=None)
equivalence_ratios = [0.75, 1.0, 1.5]
h2_fractions = np.linspace(0.1, 0.5, 5)

results = {phi: {'H2%': [], 'Heat_kJ': [], 'NOx_ppm': [], 'Final_Temp': [],
                 'TimeHistories': {}} for phi in equivalence_ratios}

# Run simulations
for phi in equivalence_ratios:
    print(f"\n--- φ = {phi} ---")
    for h2_frac in h2_fractions:
        c3h8_frac = 1 - h2_frac  # Propane fraction
        fuel_mixture = f'C3H8:{c3h8_frac}, H2:{h2_frac}'
        gas.TP = 300.0, ct.one_atm
        gas.set_equivalence_ratio(phi, fuel=fuel_mixture, oxidizer='O2:1.0, N2:3.76')

        inlet = ct.Reservoir(gas)
        gas.equilibrate('HP')
        combustor = ct.IdealGasReactor(gas)
        combustor.volume = 1.0
        exhaust = ct.Reservoir(gas)

        def mdot(t): return combustor.mass / residence_time
        inlet_mfc = ct.MassFlowController(inlet, combustor, mdot=mdot)
        outlet_mfc = ct.PressureController(combustor, exhaust, primary=inlet_mfc, K=0.01)

        sim = ct.ReactorNet([combustor])
        residence_time = 0.1

        times, heat_rates, nox_vals, temps = [], [], [], []

        while combustor.T > 500:
            sim.initial_time = 0.0
            sim.advance_to_steady_state()
            times.append(residence_time)
            heat_rates.append(combustor.thermo.heat_release_rate)
            temps.append(combustor.T)
            NOx_ppm = (combustor.thermo['NO'].X[0] + combustor.thermo['NO2'].X[0]) * 1e6
            nox_vals.append(NOx_ppm)
            residence_time *= 0.9

        Q_kJ = np.trapz(heat_rates, times) * combustor.volume / 1000
        NOx_peak = max(nox_vals)
        T_final = combustor.T
        H2_percent = round(h2_frac * 100, 1)

        results[phi]['H2%'].append(H2_percent)
        results[phi]['Heat_kJ'].append(Q_kJ)
        results[phi]['NOx_ppm'].append(NOx_peak)
        results[phi]['Final_Temp'].append(T_final)
        results[phi]['TimeHistories'][H2_percent] = {
            'time': times, 'temp': temps, 'heat': heat_rates
        }

        print(f"H₂ {H2_percent:.0f}% → Heat: {Q_kJ:.2f} kJ, Temp: {T_final:.1f} K, NOx: {NOx_peak:.1f} ppm")

# === Original Summary Plots ===

# Total Heat Released vs H₂%
plt.figure(figsize=(8, 5))
for phi in equivalence_ratios:
    plt.plot(results[phi]['H2%'], results[phi]['Heat_kJ'], marker='o', label=f'ϕ = {phi}')
plt.xlabel('Hydrogen %')
plt.ylabel('Total Heat Released (kJ)')
plt.title('Total Heat Released vs Hydrogen % (All φ)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# NOx vs H₂%
plt.figure(figsize=(8, 5))
for phi in equivalence_ratios:
    plt.plot(results[phi]['H2%'], results[phi]['NOx_ppm'], marker='s', label=f'ϕ = {phi}')
plt.xlabel('Hydrogen %')
plt.ylabel('NOx Emissions (ppm)')
plt.title('NOx Emissions vs Hydrogen % (All φ)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Final Temperature vs H₂%
plt.figure(figsize=(8, 5))
for phi in equivalence_ratios:
    plt.plot(results[phi]['H2%'], results[phi]['Final_Temp'], marker='^', label=f'ϕ = {phi}')
plt.xlabel('Hydrogen %')
plt.ylabel('Final Temperature (K)')
plt.title('Final Temperature vs Hydrogen % (All φ)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Heat vs Final Temperature
plt.figure(figsize=(8, 5))
for phi in equivalence_ratios:
    plt.plot(results[phi]['Final_Temp'], results[phi]['Heat_kJ'], marker='d', label=f'ϕ = {phi}')
plt.xlabel('Final Temperature (K)')
plt.ylabel('Total Heat Released (kJ)')
plt.title('Heat Released vs Final Temperature (per φ)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# === Temperature Rise and Heat Release Rate vs Time ===

for phi in equivalence_ratios:
    plt.figure(figsize=(8, 5))
    for H2_percent, data in results[phi]['TimeHistories'].items():
        plt.plot(data['time'], data['temp'], marker='o', label=f"H₂ {H2_percent}%")
    plt.xscale('log')
    plt.xlabel("Residence Time [s] (log scale)")
    plt.ylabel("Temperature [K]")
    plt.title(f"Temperature Rise vs Time (ϕ = {phi})")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(8, 5))
    for H2_percent, data in results[phi]['TimeHistories'].items():
        plt.plot(data['time'], data['heat'], marker='s', label=f"H₂ {H2_percent}%")
    plt.xscale('log')
    plt.xlabel("Residence Time [s] (log scale)")
    plt.ylabel("Heat Release Rate [W/m³]")
    plt.title(f"Heat Release Rate vs Time (ϕ = {phi})")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


# In[2]:


import numpy as np
import matplotlib.pyplot as plt

# Input data
hydrogen_percent = np.array([10, 20, 30, 40, 50])
phi_values = np.array([0.75, 1.0, 1.5])

# Manually extracted NOx data from propane-hydrogen chart
nox_data = np.array([
    [110, 115, 122, 130, 145],     # φ = 0.75
    [860, 900, 950, 1010, 1100],   # φ = 1.0
    [80, 82, 85, 88, 90]           # φ = 1.5
])

# Manually extracted Heat Release data from propane-hydrogen chart
heat_data = np.abs(np.array([
    [2470, 2520, 2570, 2650, 2720],    # φ = 0.75
    [2820, 2870, 2900, 2950, 3050],    # φ = 1.0
    [2530, 2570, 2610, 2700, 2750]     # φ = 1.5
]))

# Create meshgrid
X, Y = np.meshgrid(hydrogen_percent, phi_values)

# ---- Heat Map for NOx ----
plt.figure(figsize=(8, 5))
nox_plot = plt.contourf(X, Y, nox_data, levels=30, cmap='plasma')
plt.colorbar(nox_plot, label='NOx (ppm)')
plt.title('Heat Map of NOx Emissions (Propane + Hydrogen)')
plt.xlabel('Hydrogen %')
plt.ylabel('Equivalence Ratio (ϕ)')
plt.grid(True, linestyle='--', alpha=0.3)
plt.tight_layout()
plt.show()

# ---- Heat Map for Heat Release ----
plt.figure(figsize=(8, 5))
heat_plot = plt.contourf(X, Y, heat_data, levels=30, cmap='inferno')
plt.colorbar(heat_plot, label='Heat Released (kJ)')
plt.title('Heat Map of Total Heat Released (Propane + Hydrogen)')
plt.xlabel('Hydrogen %')
plt.ylabel('Equivalence Ratio (ϕ)')
plt.grid(True, linestyle='--', alpha=0.3)
plt.tight_layout()
plt.show()


# In[3]:


import numpy as np
import matplotlib.pyplot as plt
import cantera as ct

# Fuels and equivalence ratios
fuels = ['CH4', 'H2']
phis = [0.75, 1.0, 1.5]

# Data storage
data = {fuel: {phi: {'tres': [], 'T': [], 'Q': [], 'NOx': []} for phi in phis} for fuel in fuels}
avg_data = {fuel: {'avg_Q_kJ': [], 'avg_NOx': []} for fuel in fuels}

# Simulation loop
for fuel in fuels:
    for phi in phis:
        gas = ct.Solution('gri30.yaml', transport_model=None)
        gas.TP = 300.0, ct.one_atm
        gas.set_equivalence_ratio(phi, f'{fuel}:1.0', 'O2:1.0, N2:3.76')
        
        inlet = ct.Reservoir(gas)
        gas.equilibrate('HP')
        combustor = ct.IdealGasReactor(gas)
        combustor.volume = 1.0
        exhaust = ct.Reservoir(gas)

        def mdot(t): return combustor.mass / residence_time
        inlet_mfc = ct.MassFlowController(inlet, combustor, mdot=mdot)
        outlet_mfc = ct.PressureController(combustor, exhaust, primary=inlet_mfc, K=0.01)

        sim = ct.ReactorNet([combustor])
        residence_time = 0.01

        while combustor.T > 500:
            sim.initial_time = 0.0
            sim.advance_to_steady_state()
            data[fuel][phi]['tres'].append(residence_time)
            data[fuel][phi]['T'].append(combustor.T)
            data[fuel][phi]['Q'].append(combustor.thermo.heat_release_rate)
            NOx_ppm = (combustor.thermo['NO'].X[0] + combustor.thermo['NO2'].X[0]) * 1e6
            data[fuel][phi]['NOx'].append(NOx_ppm)
            residence_time *= 0.9

        # Compute average heat (in kJ): Q_avg * tres_avg * volume / 1000
        avg_q = np.mean(data[fuel][phi]['Q'])
        avg_tres = np.mean(data[fuel][phi]['tres'])
        Q_kJ = avg_q * avg_tres * 1.0 / 1000
        avg_data[fuel]['avg_Q_kJ'].append(Q_kJ)
        avg_data[fuel]['avg_NOx'].append(np.mean(data[fuel][phi]['NOx']))

# === Plot Temperature Rise ===
for fuel in fuels:
    plt.figure(figsize=(8, 5))
    for phi in phis:
        plt.plot(data[fuel][phi]['tres'], data[fuel][phi]['T'], marker='o', label=f'ϕ = {phi}')
    plt.xscale('log')
    plt.xlabel("Residence Time [s]")
    plt.ylabel("Temperature [K]")
    plt.title(f"{fuel} - Temperature Rise vs Residence Time")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

# === Plot Heat Release Rate ===
for fuel in fuels:
    plt.figure(figsize=(8, 5))
    for phi in phis:
        plt.plot(data[fuel][phi]['tres'], data[fuel][phi]['Q'], marker='s', label=f'ϕ = {phi}')
    plt.xscale('log')
    plt.xlabel("Residence Time [s]")
    plt.ylabel("Heat Release Rate [W/m³]")
    plt.title(f"{fuel} - Heat Release Rate vs Residence Time")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

# === Plot Average Heat Release (kJ) Bar Chart ===
fig, ax = plt.subplots(figsize=(8, 5))
width = 0.35
x = np.arange(len(phis))
for i, fuel in enumerate(fuels):
    ax.bar(x + i * width, avg_data[fuel]['avg_Q_kJ'], width, label=fuel)
ax.set_xticks(x + width / 2)
ax.set_xticklabels([f'ϕ = {phi}' for phi in phis])
ax.set_ylabel('Avg. Heat Released [kJ]')
ax.set_title('Average Heat Released vs Equivalence Ratio')
ax.legend()
plt.tight_layout()
plt.show()

# === Plot Average NOx Bar Chart ===
fig, ax = plt.subplots(figsize=(8, 5))
for i, fuel in enumerate(fuels):
    ax.bar(x + i * width, avg_data[fuel]['avg_NOx'], width, label=fuel)
ax.set_xticks(x + width / 2)
ax.set_xticklabels([f'ϕ = {phi}' for phi in phis])
ax.set_ylabel('Avg. NOx [ppm]')
ax.set_title('Average NOx Emissions vs Equivalence Ratio')
ax.legend()
plt.tight_layout()
plt.show()


# In[4]:


import numpy as np
import matplotlib.pyplot as plt

# Fuel properties
LHV_fuels = {
    'CH4': 50.0,     # MJ/kg
    'CH3OH': 19.9,   # MJ/kg
    'C3H8': 46.3     # MJ/kg
}

MW_fuels = {
    'CH4': 16.04,
    'CH3OH': 32.04,
    'C3H8': 44.1
}

# Hydrogen properties
LHV_H2 = 120.0  # MJ/kg
MW_H2 = 2.016
AFR_H2 = 34.3

# Stoichiometric AFRs for base fuels (kg air/kg fuel)
AFR_fuels = {
    'CH4': 17.2,
    'CH3OH': 6.47,
    'C3H8': 15.6
}

# Hydrogen mole fractions to blend
h2_mole_fracs = np.linspace(0.1, 0.5, 5)
equiv_ratios = [0.75, 1.0, 1.5]
fuels = ['CH4', 'CH3OH', 'C3H8']

# Storage
results = {fuel: {phi: {'H2%': [], 'LHV_mix': [], 'Q_air': []} for phi in equiv_ratios} for fuel in fuels}

# Loop over fuels, phi, and H2 %
for fuel in fuels:
    for phi in equiv_ratios:
        for h2_x in h2_mole_fracs:
            fuel_x = 1.0 - h2_x

            # Convert mole fractions to mass fractions
            mass_H2 = h2_x * MW_H2
            mass_fuel = fuel_x * MW_fuels[fuel]
            total_mass = mass_H2 + mass_fuel
            Y_H2 = mass_H2 / total_mass
            Y_fuel = mass_fuel / total_mass

            # Weighted stoichiometric AFR
            afr_stoich = 1.0 / ((Y_H2 / AFR_H2) + (Y_fuel / AFR_fuels[fuel]))

            # Actual AFR
            afr_actual = afr_stoich / phi

            # Mixture LHV (MJ/kg fuel)
            LHV_mix = Y_H2 * LHV_H2 + Y_fuel * LHV_fuels[fuel]

            # Heat per kg air (kJ/kg air)
            Q_air = (LHV_mix / afr_actual) * 1000  # Convert MJ → kJ

            # Store
            results[fuel][phi]['H2%'].append(h2_x * 100)
            results[fuel][phi]['LHV_mix'].append(LHV_mix)
            results[fuel][phi]['Q_air'].append(Q_air)

# === PLOTS ===

for fuel in fuels:
    # Plot 1: LHV of mixture (MJ/kg fuel)
    plt.figure(figsize=(8, 5))
    for phi in equiv_ratios:
        plt.plot(results[fuel][phi]['H2%'], results[fuel][phi]['LHV_mix'], marker='o', label=f'ϕ = {phi}')
    plt.xlabel('H₂ % in Blend')
    plt.ylabel('LHV of Mixture [MJ/kg fuel]')
    plt.title(f'{fuel} + H₂ - LHV of Fuel Mixture')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Plot 2: Heat Release per kg air
    plt.figure(figsize=(8, 5))
    for phi in equiv_ratios:
        plt.plot(results[fuel][phi]['H2%'], results[fuel][phi]['Q_air'], marker='s', label=f'ϕ = {phi}')
    plt.xlabel('H₂ % in Blend')
    plt.ylabel('Heat Release [kJ/kg air]')
    plt.title(f'{fuel} + H₂ - Heat Release per kg Air')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


# In[5]:


import cantera as ct
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Parameters
initial_T = 300.0
initial_P = ct.one_atm
phi = 1.0
h2_fracs = np.linspace(0.1, 0.5, 5)
spark_temps = [1100.0, 1150.0, 1200.0]
residence_times = [0.001, 0.003, 0.005]
mechanism = 'gri30.yaml'

data = []

# Functions
def calculate_nox(gas):
    return sum([
        gas.X[gas.species_index(X)] if X in gas.species_names else 0.0
        for X in ['NO', 'NO2', 'N2O']
    ]) * 1e6

def get_emissions(gas):
    CO = gas.X[gas.species_index('CO')] * 1e6 if 'CO' in gas.species_names else 0.0
    CO2 = gas.X[gas.species_index('CO2')] * 1e6 if 'CO2' in gas.species_names else 0.0
    return CO, CO2

# Main simulation loop
for h2_frac in h2_fracs:
    ch4_frac = 1 - h2_frac
    fuel = f'H2:{h2_frac}, CH4:{ch4_frac}'
    oxidizer = 'O2:1.0, N2:3.76'

    # Adiabatic Equilibrium
    gas_eq = ct.Solution(mechanism)
    gas_eq.set_equivalence_ratio(phi, fuel=fuel, oxidizer=oxidizer)
    gas_eq.TP = initial_T, initial_P
    gas_eq.equilibrate('HP')
    T_ad = gas_eq.T
    NOx_ad = calculate_nox(gas_eq)
    CO_ad, CO2_ad = get_emissions(gas_eq)

    data.append({
        'H2%': h2_frac * 100,
        'Spark_Temp (K)': 'Equilibrium',
        'Residence_Time (s)': 'Equilibrium',
        'Final_Temp (K)': T_ad,
        'Delta_T (K)': T_ad - initial_T,
        'NOx (ppm)': NOx_ad,
        'CO (ppm)': CO_ad,
        'CO2 (ppm)': CO2_ad
    })

    # Simulated spark + time cases
    for T_spark in spark_temps:
        for tau in residence_times:
            gas = ct.Solution(mechanism)
            gas.set_equivalence_ratio(phi, fuel=fuel, oxidizer=oxidizer)
            gas.TP = T_spark, initial_P

            reactor = ct.IdealGasReactor(gas)
            sim = ct.ReactorNet([reactor])
            t = 0.0
            while t < tau:
                t = sim.step()

            T_final = reactor.T
            delta_T = T_final - initial_T
            NOx_ppm = calculate_nox(reactor.thermo)
            CO_ppm, CO2_ppm = get_emissions(reactor.thermo)

            data.append({
                'H2%': h2_frac * 100,
                'Spark_Temp (K)': T_spark,
                'Residence_Time (s)': tau,
                'Final_Temp (K)': T_final,
                'Delta_T (K)': delta_T,
                'NOx (ppm)': NOx_ppm,
                'CO (ppm)': CO_ppm,
                'CO2 (ppm)': CO2_ppm
            })

# Create DataFrame
df = pd.DataFrame(data)

# Split spark & equilibrium
df_eq = df[df['Spark_Temp (K)'] == 'Equilibrium']
df_spark = df[df['Spark_Temp (K)'] != 'Equilibrium']
df_spark['Spark_Temp (K)'] = df_spark['Spark_Temp (K)'].astype(float)

# === Plot NOx for each residence time ===
for tau in sorted(df_spark['Residence_Time (s)'].unique()):
    plt.figure(figsize=(8, 5))
    for T_spark in sorted(df_spark['Spark_Temp (K)'].unique()):
        subset = df_spark[
            (df_spark['Spark_Temp (K)'] == T_spark) &
            (df_spark['Residence_Time (s)'] == tau)
        ]
        label = f"{int(T_spark)}K"
        plt.plot(subset['H2%'], subset['NOx (ppm)'], marker='o', label=label)
    
    # Add equilibrium line
    plt.plot(df_eq['H2%'], df_eq['NOx (ppm)'], 'k--', label='Equilibrium', linewidth=2)

    plt.title(f"NOx vs H₂% @ {int(tau*1000)} ms")
    plt.xlabel("Hydrogen in Fuel (%)")
    plt.ylabel("NOx Emissions (ppm)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

# === Plot Final Temperature for each residence time ===
for tau in sorted(df_spark['Residence_Time (s)'].unique()):
    plt.figure(figsize=(8, 5))
    for T_spark in sorted(df_spark['Spark_Temp (K)'].unique()):
        subset = df_spark[
            (df_spark['Spark_Temp (K)'] == T_spark) &
            (df_spark['Residence_Time (s)'] == tau)
        ]
        label = f"{int(T_spark)}K"
        plt.plot(subset['H2%'], subset['Final_Temp (K)'], marker='s', label=label)

    # Add equilibrium line
    plt.plot(df_eq['H2%'], df_eq['Final_Temp (K)'], 'k--', label='Equilibrium', linewidth=2)

    plt.title(f"Final Temperature vs H₂% @ {int(tau*1000)} ms")
    plt.xlabel("Hydrogen in Fuel (%)")
    plt.ylabel("Final Temperature (K)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


# In[7]:


import matplotlib.pyplot as plt

# Hydrogen percentages in fuel (mole %)
h2_percentages = [10, 20, 30, 40, 50]

# Corresponding heat released values in MJ/m³ (example values matching the plot)
heat_ch3oh = [3.61, 3.60, 3.59, 3.58, 3.57]
heat_ch4   = [3.43, 3.43, 3.43, 3.43, 3.43]
heat_c3h8  = [3.35, 3.35, 3.35, 3.35, 3.36]

# Plotting
plt.figure(figsize=(8, 5))
plt.plot(h2_percentages, heat_ch3oh, marker='o', color='gold', label='H₂ + CH₃OH')
plt.plot(h2_percentages, heat_ch4, marker='s', color='orange', label='H₂ + CH₄')
plt.plot(h2_percentages, heat_c3h8, marker='^', color='crimson', label='H₂ + C₃H₈')

plt.title('Heat Released per 1 m³ vs. Hydrogen % in Fuel Mixture-Thermodynamic Model')
plt.xlabel('Hydrogen in Fuel (%)')
plt.ylabel('Heat Released (MJ/m³)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()


# In[ ]:




