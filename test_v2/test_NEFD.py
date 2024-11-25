# test_detector_nefd.py

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
nefd_results = {}

# Loop through each detector type to compute NEFD
for detector_type in detector_types:
    # Use the factory to create the appropriate detector instance
    detector = DetectorFactory.create_detector(detector_type)

    # Calculate NEFD using the detector's characteristics
    nefd = detector.calculate_nefd(
        aperture_size=aperture_size,
        integration_time=integration_time,
        temperature=temperature,
        background_flux=background_flux
    )

    # Save the results
    nefd_results[detector_type] = nefd

# Save NEFD results to a JSON file for further analysis
results_file_path_nefd = results_file_path / "detector_nefd_results.json"

with open(results_file_path_nefd, 'w') as results_file:
    json.dump(nefd_results, results_file, indent=4)

print(f"Results saved to {results_file_path_nefd}")

# Plot the NEFD results for visual validation
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))
plt.bar(nefd_results.keys(), nefd_results.values(), color=["#0072B2", "#E69F00", "#56B4E9"])
plt.xlabel("Detector Type")
plt.ylabel("NEFD (W/m^2/sr)")
plt.title("Noise Equivalent Flux Density (NEFD) for Different Detector Types")
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()
