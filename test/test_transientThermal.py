# example_transient_thermal.py

from libraries.physicslib import TransientThermalModel, EquilibriumTemperature
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def main():
    # Balloon parameters
    diameter = 3.0  # meters
    radius = diameter / 2
    A_S = 4 * np.pi * radius ** 2       # Surface area [m^2]
    A_C = np.pi * radius ** 2           # Cross-sectional area [m^2]

    # Heat capacity for thin balloons (~0.5 kg)
    mass_thin = 0.5  # kg
    c_p_aluminum = 904    # J/kg·K
    c_p_mylar = 1150      # J/kg·K
    mass_aluminum = mass_thin * (0.19 / (0.19 + 0.25))
    mass_mylar = mass_thin * (0.25 / (0.19 + 0.25))
    C_thin = mass_aluminum * c_p_aluminum + mass_mylar * c_p_mylar  # [J/K]

    # Heat capacity for thick balloons (~3 kg)
    mass_thick = 3.0  # kg
    mass_aluminum_thick = mass_thick * 0.5  # Assuming equal mass of aluminum and mylar
    mass_mylar_thick = mass_thick * 0.5
    C_thick = mass_aluminum_thick * c_p_aluminum + mass_mylar_thick * c_p_mylar  # [J/K]

    # Common parameters
    params_common = {
        'surface_area': A_S,
        'cross_sectional_area': A_C,
        'solar_flux': 1360.0,  # W/m^2
        'albedo': 0.3,
        'earth_ir_flux': 240.0,  # W/m^2
    }

    # Initial temperature
    T0 = 300.0  # K

    # Time span for simulation (0 to 600 seconds)
    t_span = (0, 600)
    t_eval = np.linspace(*t_span, 600)  # 1-second intervals

    # List of balloons with different coatings
    balloons = [
        {
            'label': 'Shiny Aluminum Foil',
            'alpha_v': 0.192,
            'epsilon_ir': 0.036,
            'equilibrium_temp': 454,
            'heat_capacity': C_thin,
            'mass': mass_thin
        },
        {
            'label': 'Black Paint',
            'alpha_v': 0.975,
            'epsilon_ir': 0.874,
            'equilibrium_temp': 314,
            'heat_capacity': C_thin,
            'mass': mass_thin
        },
        {
            'label': 'Aluminum Silicone Paint',
            'alpha_v': 0.25,
            'epsilon_ir': 0.28,
            'equilibrium_temp': 299,
            'heat_capacity': C_thin,
            'mass': mass_thin
        },
        {
            'label': 'White TiO2 Paint',
            'alpha_v': 0.19,
            'epsilon_ir': 0.94,
            'equilibrium_temp': 227,
            'heat_capacity': C_thin,
            'mass': mass_thin
        },


    ]

    # List of balloons with different coatings
    balloons_thick = [
        {
            'label': 'Shiny Aluminum Foil',
            'alpha_v': 0.192,
            'epsilon_ir': 0.036,
            'equilibrium_temp': 454,
            'heat_capacity': C_thick,
            'mass': mass_thick
        },
        {
            'label': 'Black Paint',
            'alpha_v': 0.975,
            'epsilon_ir': 0.874,
            'equilibrium_temp': 314,
            'heat_capacity': C_thick,
            'mass': mass_thick
        },
        {
            'label': 'Aluminum Silicone Paint',
            'alpha_v': 0.25,
            'epsilon_ir': 0.28,
            'equilibrium_temp': 299,
            'heat_capacity': C_thick,
            'mass': mass_thick
        },
        {
            'label': 'White TiO2 Paint',
            'alpha_v': 0.19,
            'epsilon_ir': 0.94,
            'equilibrium_temp': 227,
            'heat_capacity': C_thick,
            'mass': mass_thick
        },

    ]
    # Simulate temperature over time for each balloon
    for balloon in balloons_thick:
        # Parameters specific to the balloon
        params = params_common.copy()
        params.update({
            'alpha_v': balloon['alpha_v'],
            'epsilon_ir': balloon['epsilon_ir'],
            'heat_capacity': balloon['heat_capacity']
        })

        # Solve the ODE
        sol = solve_ivp(
            fun=lambda t, T: TransientThermalModel.temperature_change(T, t, params),
            t_span=t_span,
            y0=[T0],
            t_eval=t_eval,
            method='RK45'
        )

        # Extract results
        T_t = sol.y[0]
        t = sol.t

        # Print the time to reach equilibrium temperature within a tolerance
        T_eq = balloon['equilibrium_temp']
        idx_eq = np.where(np.abs(T_t - T_eq) <= 1.0)[0]
        if len(idx_eq) > 0:
            time_to_eq = t[idx_eq[0]]
        else:
            time_to_eq = t[-1]
        print(f"Balloon '{balloon['label']}' reaches equilibrium temperature ({T_eq} K) in approximately {time_to_eq:.1f} seconds.")

        # Plot the temperature over time (optional, but code is provided for future use)
        plt.plot(t, T_t, label=f"{balloon['label']} ({balloon['mass']} kg)")

    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')
    plt.title('Temperature vs. Time for Balloons with Different Coatings')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()
