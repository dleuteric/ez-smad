# example_equilibrium_temperature.py

from libraries.physicslib import EquilibriumTemperature
import pandas as pd

def main():
    # Material properties from Table A-1
    # Material properties from Table A-1
    materials = [
        {
            "Coating": "White TiO2 Paint",
            "alpha_v": 0.19,
            "epsilon_ir": 0.94,
            "alpha_epsilon_ratio": 0.202,  # 0.19 / 0.94
            "expected_temp": 227
        },
        {
            "Coating": "White Epoxy Paint",
            "alpha_v": 0.248,
            "epsilon_ir": 0.924,
            "alpha_epsilon_ratio": 0.268,  # 0.248 / 0.924
            "expected_temp": 237
        },
        {
            "Coating": "White Enamel Paint",
            "alpha_v": 0.252,
            "epsilon_ir": 0.853,
            "alpha_epsilon_ratio": 0.295,  # 0.252 / 0.853
            "expected_temp": 241
        },
        {
            "Coating": "Mylar",
            "alpha_v": 0.171,
            "epsilon_ir": 0.5,
            "alpha_epsilon_ratio": 0.342,  # 0.171 / 0.5
            "expected_temp": 247
        },
        {
            "Coating": "Aluminum Silicone Paint",
            "alpha_v": 0.25,
            "epsilon_ir": 0.28,
            "alpha_epsilon_ratio": 0.893,  # 0.25 / 0.28
            "expected_temp": 299
        },
        {
            "Coating": "Gray TiO2 Paint",
            "alpha_v": 0.87,
            "epsilon_ir": 0.87,
            "alpha_epsilon_ratio": 1.0,  # 0.87 / 0.87
            "expected_temp": 307
        },
        {
            "Coating": "Black Paint",
            "alpha_v": 0.975,
            "epsilon_ir": 0.874,
            "alpha_epsilon_ratio": 1.116,  # 0.975 / 0.874
            "expected_temp": 314
        },
        {
            "Coating": "Aluminum Paint",
            "alpha_v": 0.54,
            "epsilon_ir": 0.45,
            "alpha_epsilon_ratio": 1.2,  # 0.54 / 0.45
            "expected_temp": 320
        },
        {
            "Coating": "Aquadag Paint",
            "alpha_v": 0.782,
            "epsilon_ir": 0.49,
            "alpha_epsilon_ratio": 1.596,  # 0.782 / 0.49
            "expected_temp": 341
        },
        {
            "Coating": "Shiny Aluminum Foil",
            "alpha_v": 0.192,
            "epsilon_ir": 0.036,
            "alpha_epsilon_ratio": 5.333,  # 0.192 / 0.036
            "expected_temp": 454
        },
        {
            "Coating": "Polished Gold Plate",
            "alpha_v": 0.301,
            "epsilon_ir": 0.028,
            "alpha_epsilon_ratio": 10.75,  # 0.301 / 0.028
            "expected_temp": 540
        }
    ]

    # List to store results
    results = []

    # Loop through each material and compute equilibrium temperature
    for material in materials:
        coating = material["Coating"]
        alpha_v = material["alpha_v"]
        epsilon_ir = material["epsilon_ir"]
        expected_temp = material["expected_temp"]

        # Calculate equilibrium temperature
        T_eq = EquilibriumTemperature.equilibrium_temperature(
            alpha_v=alpha_v,
            epsilon_ir=epsilon_ir
        )

        # Append results
        results.append({
            "Coating": coating,
            "alpha_v": alpha_v,
            "epsilon_ir": epsilon_ir,
            "Expected Temp (K)": expected_temp,
            "Calculated Temp (K)": round(T_eq, 1),
            "Difference (K)": round(T_eq - expected_temp, 1)
        })

    # Create a DataFrame to display the results
    df = pd.DataFrame(results)
    print(df)

# example_equilibrium_temperature.py (add the following code)

    # Calculate equilibrium temperature in Earth's shadow
    T_eq_shadow = EquilibriumTemperature.equilibrium_temperature_in_shadow(
        epsilon_ir=1.0  # ε_IR cancels out
    )
   # print(f"\nEquilibrium Temperature in Earth's Shadow: {T_eq_shadow:.1f} K")


if __name__ == "__main__":
    main()
