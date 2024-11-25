# debug_detection_range_with_logging.py

import numpy as np
from libraries.physicslib import BandRadiation, DetectionPhysics
from libraries.detectorlib import DetectorFactory

# Test Inputs for Compatibility Testing
temperatures = np.linspace(180, 400, 5)  # Temperature values ranging from 180 K to 400 K
aperture_size = 0.5  # Aperture size in meters (from benchmark case)
integration_time = 0.01  # Integration time in seconds
cross_section = 4  # Fixed cross-sectional area in m^2
snr_values = [5, 10]  # Signal-to-noise ratios for testing
background_flux = 1e-15  # Background flux for NEFD calculation (arbitrary value)

# Create detectors using the DetectorFactory with parameters from datasheets
detectors = {
    "SWIR": DetectorFactory.create_detector("SWIR"),
    "MWIR": DetectorFactory.create_detector("MWIR"),
    "LWIR": DetectorFactory.create_detector("LWIR")
}

# Define wavelength ranges for radiant exitance calculations (updated)
# These wavelength bands are defined by the detector datasheets
wavelength_bands = {
    "SWIR": (detectors["SWIR"].wavelength_min, detectors["SWIR"].wavelength_max),
    "MWIR": (detectors["MWIR"].wavelength_min, detectors["MWIR"].wavelength_max),
    "LWIR": (detectors["LWIR"].wavelength_min, detectors["LWIR"].wavelength_max)
}

# Debug Information
print("\nDebugging Information:\n")

# Store results for plotting
radiant_exitance_results = {band_name: [] for band_name in detectors}
detection_range_results = {band_name: {snr: [] for snr in snr_values} for band_name in detectors}

# Calculate radiant exitance, NEFD, and detection range for each temperature and band
for temp in temperatures:
    print(f"Temperature: {temp} K")
    for band_name, detector in detectors.items():
        wavelength_min, wavelength_max = wavelength_bands[band_name]

        # Step 1: Calculate band radiant exitance using updated wavelength limits
        radiant_exitance = BandRadiation.band_radiant_exitance(
            temp, wavelength_min, wavelength_max
        )
        radiant_exitance_results[band_name].append(radiant_exitance)
        print(f"{band_name} Band -> Radiant Exitance: {radiant_exitance:.6e} W/m^2")

        # Step 2: Calculate NEFD using the updated detector model and datasheet values
        nefd = detector.calculate_nefd(aperture_size, integration_time, temp, background_flux)
        print(f"{band_name} Band -> NEFD: {nefd:.6e} W/m^2")

        # Step 3: Calculate power emitted using radiant exitance and cross-sectional area
        power_emitted = radiant_exitance * cross_section
        print(f"{band_name} Band -> Power Emitted: {power_emitted:.6e} W")

        # Step 4: Calculate detection range for each SNR value using updated detector information
        for snr in snr_values:
            f_number = detector.f_number
            netd = detector.netd

            try:
                detection_range = DetectionPhysics.detection_range(
                    power_emitted=power_emitted,
                    nefd_w_m2=nefd,
                    snr=snr,
                    band_radiant_exitance=radiant_exitance,
                    f_number=f_number,
                    netd=netd
                )
                print(f"{band_name} Band -> Detection Range for SNR {snr}: {detection_range / 1000:.2f} km")
                detection_range_results[band_name][snr].append(detection_range / 1000)  # Convert to kilometers
            except Exception as e:
                print(f"Error in detection range calculation for {band_name} detector at {temp} K: {e}")
                detection_range_results[band_name][snr].append(None)

# Summary Output
print("\nSummary Debug Information: Radiant Exitance, NEFD, and Detection Range Calculations")
for band_name in detectors:
    print(f"\n{band_name} Band:")
    for temp, radiant_exitance in zip(temperatures, radiant_exitance_results[band_name]):
        print(f"Temperature: {temp} K -> Radiant Exitance: {radiant_exitance:.6e} W/m^2")
    for snr in snr_values:
        print(f"  Detection Range for SNR {snr}: {detection_range_results[band_name][snr]}")
