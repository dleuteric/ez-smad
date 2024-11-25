# test_radiant_exitance.py

import numpy as np
from libraries.physicslib import BandRadiation

# Test Inputs
temperatures = [180.0, 235.0, 290.0, 345.0, 400.0]  # Temperature values in Kelvin
wavelength_ranges = {
    "SWIR": (1e-6, 3e-6),   # Short-Wave Infrared: 1 to 3 micrometers
    "MWIR": (3e-6, 8e-6),   # Medium-Wave Infrared: 3 to 8 micrometers
    "LWIR": (8e-6, 15e-6)   # Long-Wave Infrared: 8 to 15 micrometers
}
path_length = 10.0  # Path length in kilometers
water_vapor_content = 7.5  # Water vapor content in g/m^3

# Run calculations for each temperature and wavelength band
for band_name, (wavelength_min, wavelength_max) in wavelength_ranges.items():
    print(f"\n{band_name} Band:")

    for temp in temperatures:
        # Calculate atmospheric transmittance for the wavelength range midpoint
        wavelength_mid = (wavelength_min + wavelength_max) / 2
        transmittance = BandRadiation.atmospheric_transmittance(wavelength_mid, path_length, water_vapor_content)

        # Calculate radiant exitance with atmospheric attenuation
        radiant_exitance = BandRadiation.band_radiant_exitance(temp, wavelength_min, wavelength_max, transmittance)

        # Print debug information
        print(f"  Temperature: {temp} K -> Transmittance: {transmittance:.4f} -> Radiant Exitance: {radiant_exitance:.6e} W/m^2")

