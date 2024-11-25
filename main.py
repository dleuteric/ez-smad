import numpy as np
from matplotlib import pyplot as plt

from libraries import BlackbodyRadiation, GrayBodyRadiation, DetectionPhysics, PlotLib

# Input Parameters
# Temperatures range from 180 to 300 K
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

# 1. Angular Resolution Analysis (Diffraction-Limited)
for band_name, (wavelength_min, wavelength_max) in wavelength_bands.items():
    angular_resolutions = []
    wavelength = (wavelength_min + wavelength_max) / 2  # Use the midpoint wavelength for analysis
    for aperture in aperture_sizes:
        resolution = DetectionPhysics.diffraction_limited_angular_resolution(wavelength, aperture)
        angular_resolutions.append(resolution)

    # Plot: Angular Resolution vs Aperture Size (combined in a single plot)
    PlotLib.plot_relationship(aperture_sizes, angular_resolutions, xlabel="Aperture Size (m)",
                              ylabel="Angular Resolution (radians)",
                              title="Diffraction-Limited Angular Resolution vs Aperture Size",
                              legend=f"{band_name}")

# 2. Detection Range vs Temperature Analysis (with Fixed SNR, Aperture Size)
plt.figure(figsize=(10, 6))  # Combine into a single plot for different bands
for band_name, (wavelength_min, wavelength_max) in wavelength_bands.items():
    detection_ranges_temp = []
    fixed_aperture = 0.3  # Fixed aperture size [m]
    fixed_snr = 5  # Fixed signal-to-noise ratio
    wavelength = (wavelength_min + wavelength_max) / 2  # Use the midpoint wavelength for analysis
    for temp in temperatures:
        radiant_exitance = GrayBodyRadiation.radiant_exitance(temp, emissivity)  # Total radiant exitance [W/m^2]
        power_emitted = radiant_exitance * threat_cross_sections[0]  # Total emitted power (cross-section area = 0.25 m^2)
        detection_range = DetectionPhysics.detection_range(power_emitted, nefd, fixed_snr)
        detection_ranges_temp.append(detection_range / 1000)  # Convert to kilometers

    # Add plot line for each wavelength band
    plt.plot(temperatures, detection_ranges_temp, marker='o', linestyle='-', linewidth=2, label=f"{band_name}")

plt.xlabel("Temperature (K)")
plt.ylabel("Detection Range (km)")
plt.title("Detection Range vs Temperature (Aperture = 0.3 m, Cross-Section = 0.25 m², SNR = 5)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# 3. Detection Range vs Aperture Size (Parametric Analysis wrt Temperature)
plt.figure(figsize=(10, 6))  # Combine into a single plot for different bands
for band_name, (wavelength_min, wavelength_max) in wavelength_bands.items():
    wavelength = (wavelength_min + wavelength_max) / 2  # Use the midpoint wavelength for analysis
    for idx, temp in enumerate(temperatures):
        detection_ranges_aperture = []
        for aperture in aperture_sizes:
            radiant_exitance = GrayBodyRadiation.radiant_exitance(temp, emissivity)
            power_emitted = radiant_exitance * threat_cross_sections[0]
            detection_range = DetectionPhysics.detection_range(power_emitted, nefd, fixed_snr)
            detection_ranges_aperture.append(detection_range / 1000)  # Convert to kilometers
        plt.plot(aperture_sizes, detection_ranges_aperture, marker='o', linestyle='-', linewidth=2, label=f"{band_name}, Temp = {temp} K")

plt.xlabel("Aperture Size (m)")
plt.ylabel("Detection Range (km)")
plt.title("Detection Range vs Aperture Size (Parametric wrt Temperature)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# 4. Detection Range vs Cross-Section Area (Parametric wrt SNR)
plt.figure(figsize=(10, 6))  # Combine into a single plot for different bands
for band_name, (wavelength_min, wavelength_max) in wavelength_bands.items():
    fixed_temp = 230  # Fixed temperature [K]
    wavelength = (wavelength_min + wavelength_max) / 2  # Use the midpoint wavelength for analysis
    for cross_section in threat_cross_sections:
        detection_ranges_cross_section = []
        for snr in snr_values:
            radiant_exitance = GrayBodyRadiation.radiant_exitance(fixed_temp, emissivity)
            power_emitted = radiant_exitance * cross_section
            detection_range = DetectionPhysics.detection_range(power_emitted, nefd, snr)
            detection_ranges_cross_section.append(detection_range / 1000)  # Convert to kilometers
        plt.plot(snr_values, detection_ranges_cross_section, marker='o', linestyle='-', linewidth=2, label=f"{band_name}, Cross-Section = {cross_section} m²")

plt.xlabel("Signal-to-Noise Ratio (SNR)")
plt.ylabel("Detection Range (km)")
plt.title("Detection Range vs SNR for Different Cross-Section Areas")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# 5. Comprehensive 3D Parametric Analysis
# 3D Surface Plot of Detection Range wrt Temperature and Aperture Size
for band_name, (wavelength_min, wavelength_max) in wavelength_bands.items():
    temperature_mesh, aperture_mesh = np.meshgrid(temperatures, aperture_sizes, indexing='ij')
    detection_ranges_3d = np.zeros_like(temperature_mesh)
    wavelength = (wavelength_min + wavelength_max) / 2  # Use the midpoint wavelength for analysis

    for i in range(temperature_mesh.shape[0]):
        for j in range(temperature_mesh.shape[1]):
            temp = temperature_mesh[i, j]
            aperture = aperture_mesh[i, j]
            radiant_exitance = GrayBodyRadiation.radiant_exitance(temp, emissivity)
            power_emitted = radiant_exitance * threat_cross_sections[0]  # Use fixed cross-section
            detection_ranges_3d[i, j] = DetectionPhysics.detection_range(power_emitted, nefd, fixed_snr) / 1000  # km

    # Plot the parametric analysis using 3D surface plot
    additional_params = {"Cross-Section": f"{threat_cross_sections[0]} m²", "SNR": f"{fixed_snr} (Fixed)"}
    PlotLib.plot_parametric_analysis(temperatures, aperture_sizes, detection_ranges_3d, xlabel="Temperature (K)",
                                     ylabel="Aperture Size (m)", zlabel="Detection Range (km)",
                                     additional_params=additional_params)

# Conclusion
# This main script demonstrates the utility of the tool to model detection ranges, angular resolutions, and
# temperature dependencies in a space-based infrared observation system, inspired by SBIRS-Low capabilities.
