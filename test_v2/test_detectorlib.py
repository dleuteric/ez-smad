# test_detector_characteristics.py

import numpy as np
from libraries.detectorlib import DetectorFactory
from pathlib import Path
import json

# Define the project root dynamically
PROJECT_ROOT = Path(__file__).parent

# SAVE DIRECTORIES
results_file_path = PROJECT_ROOT / "results"

# Ensure the results directory exists
results_file_path.mkdir(exist_ok=True)

# Define the detector types to be tested
detector_types = ["SWIR", "MWIR", "LWIR"]

# Define common test parameters
aperture_size = 0.3  # Aperture size in meters
integration_time = 0.01  # Integration time in seconds (e.g., 10ms)
temperature = 250  # Temperature in Kelvin
background_flux = 5e-18  # Background flux in W/m^2/sr (assumed value for test)

# Store results for validation
detector_properties = {}

# Loop through each detector type to compute properties
for detector_type in detector_types:
    # Use the factory to create the appropriate detector instance
    detector = DetectorFactory.create_detector(detector_type)

    # Save detector properties
    detector_properties[detector_type] = {
        "Quantum Efficiency": detector.get_quantum_efficiency(),
        "Dark Current Parameters": detector.dark_current_params,
        "Shot Noise Parameters": detector.shot_noise_params,
        "Thermal Noise Parameters": detector.thermal_noise_params
    }

# Save properties to a JSON file for further analysis
results_file_path_properties = results_file_path / "detector_properties_results.json"

with open(results_file_path_properties, 'w') as results_file:
    json.dump(detector_properties, results_file, indent=4)

print(f"Results saved to {results_file_path_properties}")

# Plot Quantum Efficiency for visual validation
import matplotlib.pyplot as plt

qe_results = {detector_type: props["Quantum Efficiency"] for detector_type, props in detector_properties.items()}

plt.figure(figsize=(10, 6))
plt.bar(qe_results.keys(), qe_results.values(), color=["#0072B2", "#E69F00", "#56B4E9"])
plt.xlabel("Detector Type")
plt.ylabel("Quantum Efficiency")
plt.title("Quantum Efficiency for Different Detector Types")
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()
