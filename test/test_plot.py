import numpy as np
from libraries.plotlib import PlotLib

# Set plot defaults
PlotLib.set_plot_defaults()

# 1. Detection Range vs SNR Validation (Parametric wrt Temperature)
# SNR range from 5 to 10, at different temperatures (fixed cross-section and aperture)
snr_values = np.linspace(5, 10, 6)  # SNR ranging from 5 to 10
temperatures = [180, 230, 270, 300, 330]  # Temperatures in Kelvin
detection_ranges_snr_temp = []
for temp in temperatures:
    detection_range = 1500 * np.log(snr_values) * (1 + (temp - 180) / 120)  # Creative formula for validation
    detection_ranges_snr_temp.append(detection_range)

# Plot detection range vs SNR for different temperatures
PlotLib.plot_detection_range_vs_snr_parametric(snr_values, detection_ranges_snr_temp, temperatures, fixed_parameter='temperature', parameter_value='1 m^2')

# 2. Detection Range vs SNR Validation (Parametric wrt Cross Section)
cross_sections = [0.25, 0.5, 1.0, 2.0]  # Cross section areas in m^2
detection_ranges_snr_cross_section = []
for cross_sec in cross_sections:
    detection_range = 2000 * np.log(snr_values) * (cross_sec / 1.0)  # Creative formula for validation
    detection_ranges_snr_cross_section.append(detection_range)

# Plot detection range vs SNR for different cross sections
PlotLib.plot_detection_range_vs_snr_parametric(snr_values, detection_ranges_snr_cross_section, cross_sections, fixed_parameter='cross section', parameter_value='230 K')

# 3. Aperture vs Detection Range Validation (Parametric wrt Temperature)
aperture_sizes = np.array([0.05, 0.1, 0.15, 0.2, 0.25, 0.3])  # Aperture sizes in meters
temperatures = [180, 230, 280]  # Temperatures in Kelvin
detection_ranges_aperture_temp = []
for temp in temperatures:
    detection_range = 2000 * (aperture_sizes ** 0.5) * (1 + (temp - 180) / 120)  # Creative formula for validation
    detection_ranges_aperture_temp.append(detection_range)

# Plot detection range vs aperture size for different temperatures
for idx, temp in enumerate(temperatures):
    PlotLib.plot_detection_range_vs_aperture_parametric(aperture_sizes, detection_ranges_aperture_temp[idx], fixed_param_value='SNR = 10', fixed_temp=temp)

# 4. Time Evolution of Temperature Validation
# Simulate thick balloons (mass 3 kg, radius 3 m)
# Temperature oscillation over time from initial temperature of 300 K
time_values = np.linspace(0, 120, 50)  # Time points in minutes
temperature_values = 300 + 20 * np.sin(time_values / 120 * 2 * np.pi)  # Sinusoidal temperature oscillation

# Plot time evolution of temperature for the balloon
PlotLib.plot_time_evolution_of_temperature(time_values, temperature_values, material_name="Thick Balloon Decoy", total_time_minutes=120)

# 5. Parametric Analysis Validation
# Detection range as a function of target cross-section and SNR
cross_sections = np.linspace(0.25, 5, 10)  # Cross-section areas in m^2
snr_values = np.linspace(5, 10, 5)  # SNR ranging from 5 to 10
cross_section_mesh, snr_mesh = np.meshgrid(cross_sections, snr_values)
detection_ranges_parametric = 1000 * np.sqrt(cross_section_mesh) * np.log(snr_mesh)  # Creative relationship for validation

# Plot parametric analysis with a 3D surface
additional_params = {"Fixed Aperture Size": "0.3 m", "Fixed Temperature": "230 K"}
PlotLib.plot_parametric_analysis(cross_sections, snr_values, detection_ranges_parametric, xlabel="Cross-Section Area (m^2)", ylabel="SNR", zlabel="Detection Range (km)", additional_params=additional_params)

# 6. Multiple Scenarios Validation (Parametric wrt Instrument Aperture and Target Temperature)
# Comparing detection ranges for different types of sensors (e.g., different wavebands)
temperatures = [200, 250, 300]  # Different temperatures for comparison
aperture_sizes = [0.1, 0.2, 0.3]  # Different aperture sizes for comparison
scenarios = []
for temp, aperture in zip(temperatures, aperture_sizes):
    detection_range = 1000 * (1 + (temp - 180) / 120) * (aperture / 0.1)
    scenarios.append((temp_values, detection_range, f'Temperature = {temp} K, Aperture = {aperture} m'))

# Plot multiple scenarios to compare sensor performance with parametric variation
PlotLib.plot_multiple_scenarios(scenarios, xlabel='Temperature (K)', ylabel='Detection Range (km)')
