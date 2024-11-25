# example_graybody_validation.py

from libraries.physicslib import GrayBodyRadiation, BlackbodyRadiation
from libraries.mathlib_0 import NumericalIntegration
import numpy as np

def main():
    # Given temperature and emissivity
    T = 300       # Temperature in Kelvin
    epsilon = 0.8 # Emissivity of the gray body

    # Compute total radiant exitance M for gray body
    M_gray = GrayBodyRadiation.radiant_exitance(T, epsilon)  # in W/m^2

    # Convert M_gray to W/cm^2 (since 1 W/m^2 = 1e-4 W/cm^2)
    M_gray_W_cm2 = M_gray * 1e-4

    print(f"Total radiant exitance M for gray body at {T} K and ε = {epsilon}: {M_gray_W_cm2:.6f} W/cm²")

    # Given area
    area = 4  # Area in m²

    # Compute total power emitted by gray body warhead
    total_power_gray = M_gray * area  # in W

    print(f"Total power emitted by gray body warhead: {total_power_gray:.1f} W")

    # Define wavelength range from 8 µm to 12 µm
    wavelength_min = 8e-6   # 8 µm in meters
    wavelength_max = 12e-6  # 12 µm in meters

    # Ensure the number of samples is odd for Simpson's rule
    num_samples = 1001  # Must be odd for Simpson's rule
    wavelengths = np.linspace(wavelength_min, wavelength_max, num_samples)

    # Compute spectral radiant exitance at each wavelength for gray body
    spectral_exitance_gray = np.array([
        GrayBodyRadiation.spectral_radiant_exitance(T, wl, epsilon) for wl in wavelengths
    ])  # in W/m²·m⁻¹

    # Integrate spectral exitance over the wavelength range to get M_band_gray
    M_band_gray = NumericalIntegration.simpson_rule(spectral_exitance_gray, wavelengths)  # in W/m²

    # Convert M_band_gray to W/cm²
    M_band_gray_W_cm2 = M_band_gray * 1e-4

    print(f"Radiant exitance in 8–12 µm band for gray body: {M_band_gray_W_cm2:.6f} W/cm²")

    # Compute power emitted in the band by gray body warhead
    power_band_gray = M_band_gray * area  # in W

    print(f"Power emitted in 8–12 µm band by gray body warhead: {power_band_gray:.1f} W")

    # Calculate absorptivity and reflectivity
    absorptivity = GrayBodyRadiation.absorptivity(epsilon)
    reflectivity = GrayBodyRadiation.reflectivity(epsilon)

    print(f"Absorptivity α: {absorptivity}")
    print(f"Reflectivity ρ: {reflectivity}")

if __name__ == "__main__":
    main()
