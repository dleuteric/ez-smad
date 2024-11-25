import numpy as np
from libraries import BlackbodyRadiation, GrayBodyRadiation, DetectionPhysics, PlotLib, SaveLib

# Input Parameters
temperatures = np.linspace(180, 300, 6)  # Target temperatures [K]
emissivity = 0.95  # Emissivity of the target (assumed constant)
aperture_sizes = [0.05, 0.1, 0.15, 0.2, 0.3]  # Aperture diameters in meters
wavelength_bands = {
    "SWIR": (1e-6, 3e-6),  # Short-Wave Infrared: 1 to 3 micrometers
    "MWIR": (3e-6, 6e-6),  # Medium-Wave Infrared: 3 to 6 micrometers
    "LWIR": (6e-6, 16e-6)  # Long-Wave Infrared: 6 to 16 micrometers
}
snr_values = np.linspace(1, 10, 10)  # Signal-to-noise ratios
nefd = 6e-19  # Noise Equivalent Flux Density [W/m²]
threat_cross_sections = [0.25, 1.0]  # Cross-sectional areas of the target in m^2

# Detection Range vs Target Temperature @SNR 5, varying Cross Section
# 3 curves in the same plot for SWIR, MWIR, LWIR
# Fixed SNR
fixed_snr = 500

# Store detection ranges for each wavelength band and cross-section
detection_ranges_flattened = {}

# Loop through each wavelength band
for band_name, (wavelength_min, wavelength_max) in wavelength_bands.items():
    wavelength = (wavelength_min + wavelength_max) / 2  # Use the midpoint wavelength for analysis

    # Loop through each cross-sectional area
    for cross_section in threat_cross_sections:
        detection_ranges = []

        # Calculate detection range for each temperature
        for temp in temperatures:
            radiant_exitance = GrayBodyRadiation.radiant_exitance(temp, emissivity)  # Total radiant exitance [W/m^2]
            power_emitted = radiant_exitance * cross_section  # Total emitted power (cross-section area)
            detection_range = DetectionPhysics.detection_range(power_emitted, nefd, fixed_snr)
            detection_ranges.append(detection_range / 1000)  # Convert to kilometers

        label = f"{band_name} - Cross-Section {cross_section} m²"
        detection_ranges_flattened[label] = detection_ranges

# Plot Detection Range vs Temperature for Different Wavelength Bands and Cross Sections
PlotLib.plot_multiple_relationships(temperatures, detection_ranges_flattened, xlabel="Temperature (K)",
                                    ylabel="Detection Range (km)",
                                    title=f"Detection Range vs Temperature (SNR = {fixed_snr})")

# Save Results to JSON
SaveLib.save_to_json(detection_ranges_flattened, "detection_analysis")
