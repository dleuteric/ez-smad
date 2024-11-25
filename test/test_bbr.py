# example_bbr_validation.py

from libraries.physicslib import BlackbodyRadiation as bb
from libraries.mathlib_0 import NumericalIntegration as nint
import numpy as np
from numpy import trapezoid

def main():
    # Given temperature
    T = 230  # Temperature in Kelvin

    # Compute total radiant exitance M
    M = bb.radiant_exitance(T)  # in W/m^2

    # Convert M to W/cm^2 (since 1 W/m^2 = 1e-4 W/cm^2)
    M_W_cm2 = M * 1e-4

    print(f"Total radiant exitance M at {T} K: {M_W_cm2:.6f} W/cm²")

    # Given area
    area = 0.25  # Area in m²

    # Compute total power emitted by warhead
    total_power = M * area  # in W

    print(f"Total power emitted by warhead: {total_power:.1f} W")

    # Define wavelength range from 8 µm to 12 µm
    wavelength_min = 8e-6   # 8 µm in meters
    wavelength_max = 12e-6  # 12 µm in meters

    # Generate wavelengths for integration
    wavelengths = np.linspace(wavelength_min, wavelength_max, 1001)  # 1001 points (Simpson's rule requires odd number sample)

    # Compute spectral radiant exitance at each wavelength
    spectral_exitance = np.array([
        bb.spectral_radiant_exitance(T, wl) for wl in wavelengths
    ])  # in W/m²·m⁻¹

    # Integrate spectral exitance over the wavelength range to get M_band
    M_band = nint.simpson_rule(spectral_exitance, wavelengths)  # in W/m²

    # Convert M_band to W/cm²
    M_band_W_cm2 = M_band * 1e-4

    print(f"Radiant exitance in 8–12 µm band: {M_band_W_cm2:.6f} W/cm²")

    # Compute power emitted in the band by warhead
    power_band = M_band * area  # in W

    print(f"Power emitted in 8–12 µm band by warhead: {power_band:.1f} W")

if __name__ == "__main__":
    main()
