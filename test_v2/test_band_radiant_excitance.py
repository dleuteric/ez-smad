# test_band_radiant_exitance.py

import json
import numpy as np
from libraries.physicslib import BandRadiation
from libraries.plotlib import PlotLib
from pathlib import Path

# Define the project root dynamically
PROJECT_ROOT = Path(__file__).parent

# SOURCE DIRECTORIES
LIBRARIES_DIR = PROJECT_ROOT / "libraries"

# SAVE DIRECTORIES
results_file_path = PROJECT_ROOT / "results"

# Ensure the results directory exists
results_file_path.mkdir(exist_ok=True)

# Test inputs
#temperatures = [180, 230, 270, 300, 500, 800, 1200, 1600]  # Temperatures in Kelvin for testing
temperatures = np.linspace(180, 1600, 1000)
wavelength_bands = {
    "SWIR": (1e-6, 3e-6),  # Short-Wave Infrared: 1 to 3 micrometers
    "MWIR": (3e-6, 8e-6),  # Medium-Wave Infrared: 3 to 6 micrometers
    "LWIR": (8e-6, 15e-6)  # Long-Wave Infrared: 6 to 16 micrometers
}

# Placeholder for the results
results = {}

# Calculate band radiant exitance for each temperature and wavelength band
for temp in temperatures:
    results[temp] = {}
    for band_name, (wavelength_min, wavelength_max) in wavelength_bands.items():
        radiant_exitance = BandRadiation.band_radiant_exitance(temp, wavelength_min, wavelength_max)
        results[temp][band_name] = radiant_exitance

# Save results to a JSON file for further analysis
results_file_path = Path(results_file_path / "band_radiant_exitance_results.json")

with open(results_file_path, 'w') as results_file:
    json.dump(results, results_file, indent=4)

print(f"Results saved to {results_file_path}")

# Plotting the results
PlotLib.plot_band_radiant_exitance(temperatures, results)
