# test_validation_new_range.py

import numpy as np
from libraries.detectorlib import DetectorFactory
from libraries.physicslib import DetectionPhysics, BandRadiation

# Test Inputs
snr_values = [5, 10]  # Signal-to-noise ratios for testing
cross_sections = [0.25, 1.0]  # Cross-sectional areas of the target in m^2
aperture_size = 0.3  # Aperture size in meters
integration_time = 0.1  # Integration time in seconds
temperatures = [180, 270, 300, 1600]  # Temperatures in Kelvin for testing
power_emitted = 1e-3  # Example emitted power in Watts
wavelength_bands = {
    "SWIR": (1e-6, 3e-6),  # Short-Wave Infrared: 1 to 3 micrometers
    "MWIR": (3e-6, 8e-6),  # Medium-Wave Infrared: 3 to 8 micrometers
    "LWIR": (8e-6, 15e-6)  # Long-Wave Infrared: 8 to 15 micrometers
}

# Background Flux Estimate (W/m²)
background_flux = 1e-8  # Placeholder value, should be updated based on observation scenario

# Detector Types to Test
detector_types = ["SWIR", "MWIR", "LWIR"]

# Store detection range results for validation
detection_range_results = {}

# Iterate through each detector type and temperature
for detector_type in detector_types:
    detector = DetectorFactory.create_detector(detector_type)
    detection_range_results[detector_type] = {}

    for temp in temperatures:
        detection_range_results[detector_type][f"Temperature {temp} K"] = {}

        for snr in snr_values:
            for cross_section in cross_sections:
                # Calculate emitted power using cross-section and temperature
                # Assuming power_emitted is proportional to cross_section for simplicity
                power_emitted_actual = power_emitted * cross_section

                # Calculate NEFD using the detector model
                nefd = detector.calculate_nefd(aperture_size, integration_time, temp, background_flux)

                # Calculate band radiant exitance using the BandRadiation class
                wavelength_min, wavelength_max = wavelength_bands[detector_type]
                band_radiant_exitance = BandRadiation.band_radiant_exitance(temp, wavelength_min, wavelength_max)

                # Calculate detection range using DetectionPhysics from physicslib
                detection_range = DetectionPhysics.detection_range(
                    power_emitted=power_emitted_actual,
                    nefd_w_m2=nefd,
                    snr=snr,
                    band_radiant_exitance=band_radiant_exitance
                )

                label = f"Cross-Section {cross_section} m² - SNR {snr}"
                detection_range_results[detector_type][f"Temperature {temp} K"][label] = detection_range / 1000  # Convert to kilometers

# Print detection range results for comparison
for detector_type, temp_data in detection_range_results.items():
    print(f"\nDetector Type: {detector_type}")
    for temp, ranges in temp_data.items():
        print(f"  {temp}:")
        for label, range_km in ranges.items():
            print(f"    {label}: {range_km:.2f} km")
