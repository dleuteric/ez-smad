import numpy as np
from libraries.plotlib import PlotLib
from libraries.physicslib import BlackbodyRadiation, GrayBodyRadiation, DetectionPhysics
from libraries.mathlib_0 import NumericalIntegration

# Set plot defaults
PlotLib.set_plot_defaults()

def main():
    """
    Main function to showcase the potential of the SBIRS-Low Analysis Tool using actual physics-based calculations.
    This example demonstrates the use of SBIRS-Low for detecting and tracking ballistic missile targets in space,
    utilizing realistic physics relationships.
    """
    # Constants
    SENSOR_APERTURE_DIAMETER = 0.3  # meters
    SNR_VALUES = np.linspace(5, 10, 6)  # SNR ranging from 5 to 10
    TEMPERATURES = [180, 230, 300]  # Temperatures in Kelvin
    CROSS_SECTIONS = [0.25, 0.5, 1.0, 2.0]  # Cross-section areas in m^2
    WAVELENGTH_BAND = (6e-6, 16e-6)  # Wavelength band for long-wave IR in meters
    FIXED_SNR = 10

    # 1. Detection Range vs SNR (Real Physics: Planck's Law and Emissivity)
    detection_ranges_snr_temp = []
    for temp in TEMPERATURES:
        detection_range_values = []
        for snr in SNR_VALUES:
            # Calculate emissivity and radiance
            spectral_radiance = BlackbodyRadiation.spectral_radiant_exitance(temp, WAVELENGTH_BAND[0])
            emissivity = GrayBodyRadiation.calculate_emissivity(spectral_radiance, spectral_radiance)
            # Calculate detection range based on emissivity and target properties
            effective_radiance = emissivity * spectral_radiance
            detection_range = DetectionPhysics.calculateDetectionRange(effective_radiance, sensorSensitivity=6e-19) * (snr / FIXED_SNR)
            detection_range_values.append(detection_range)
        detection_ranges_snr_temp.append(detection_range_values)

    # Plot detection range vs SNR for different temperatures
    PlotLib.plot_detection_range_vs_snr_parametric(SNR_VALUES, detection_ranges_snr_temp, TEMPERATURES, fixed_parameter='temperature', parameter_value='Cross Section = 1 m^2')

    # 2. Detection Range vs SNR (Parametric wrt Cross Section)
    detection_ranges_snr_cross_section = []
    for cross_sec in CROSS_SECTIONS:
        detection_range_values = []
        for snr in SNR_VALUES:
            # Assume a fixed temperature
            temp = 230  # Kelvin
            spectral_radiance = BlackbodyRadiation.spectralRadiantExitance(temp, WAVELENGTH_BAND[0])
            emissivity = GrayBodyRadiation.calculate_emissivity(spectral_radiance, spectral_radiance)
            effective_radiance = emissivity * spectral_radiance * cross_sec
            detection_range = DetectionPhysics.calculateDetectionRange(effective_radiance, sensorSensitivity=6e-19) * (snr / FIXED_SNR)
            detection_range_values.append(detection_range)
        detection_ranges_snr_cross_section.append(detection_range_values)

    # Plot detection range vs SNR for different cross sections
    PlotLib.plot_detection_range_vs_snr_parametric(SNR_VALUES, detection_ranges_snr_cross_section, CROSS_SECTIONS, fixed_parameter='cross section', parameter_value='Temperature = 230 K')

    # 3. Detection Range vs Aperture Size Validation
    aperture_sizes = np.array([0.05, 0.1, 0.15, 0.2, 0.25, 0.3])  # Aperture sizes in meters
    detection_ranges_aperture_temp = []
    for temp in TEMPERATURES:
        detection_range_values = []
        for aperture in aperture_sizes:
            # Use diffraction-limited resolution for sensor aperture size
            resolution = WAVELENGTH_BAND[0] / aperture
            detection_range = DetectionPhysics.calculateDetectionRange(effective_radiance, sensorSensitivity=6e-19) / resolution
            detection_range_values.append(detection_range)
        detection_ranges_aperture_temp.append(detection_range_values)

    # Plot detection range vs aperture size for different temperatures
    for idx, temp in enumerate(TEMPERATURES):
        PlotLib.plot_detection_range_vs_aperture_parametric(aperture_sizes, detection_ranges_aperture_temp[idx], fixed_param_value=f'SNR = {FIXED_SNR}', fixed_temp=temp)

    # 4. Time Evolution of Temperature Validation
    time_values = np.linspace(0, 120, 50)  # Time points in minutes
    temp_initial = 300  # Initial temperature in Kelvin
    temperature_values = [temp_initial]
    heat_capacity = 460  # J/K for balloon
    area = 4  # m^2 for the surface area
    emissivity = 0.9

    for t in time_values[1:]:
        # Calculate power emitted over time using Stefan-Boltzmann law
        power_emitted = BlackbodyRadiation.radiantEmittance(temperature_values[-1]) * area * emissivity
        temp_change = -power_emitted * (t - time_values[1]) / heat_capacity
        temperature_values.append(temperature_values[-1] + temp_change)

    # Plot time evolution of temperature for the balloon
    PlotLib.plot_time_evolution_of_temperature(time_values, temperature_values, material_name="Thick Balloon Decoy", total_time_minutes=120)

    # 5. Parametric Analysis Validation
    cross_sections = np.linspace(0.25, 5, 10)  # Cross-section areas in m^2
    snr_values = np.linspace(5, 10, 5)  # SNR ranging from 5 to 10
    cross_section_mesh, snr_mesh = np.meshgrid(cross_sections, snr_values)
    detection_ranges_parametric = np.zeros_like(cross_section_mesh)
    for i in range(len(cross_sections)):
        for j in range(len(snr_values)):
            temp = 230  # Kelvin
            spectral_radiance = BlackbodyRadiation.spectralRadiantExitance(temp, WAVELENGTH_BAND[0])
            effective_radiance = spectral_radiance * cross_sections[i] * GrayBodyRadiation.calculate_emissivity(spectral_radiance, spectral_radiance)
            detection_ranges_parametric[j, i] = DetectionPhysics.calculateDetectionRange(effective_radiance, sensorSensitivity=6e-19) * (snr_values[j] / FIXED_SNR)

    # Plot parametric analysis with a 3D surface
    additional_params = {"Fixed Aperture Size": f"{SENSOR_APERTURE_DIAMETER} m", "Fixed Temperature": "230 K"}
    PlotLib.plot_parametric_analysis(cross_sections, snr_values, detection_ranges_parametric, xlabel="Cross-Section Area (m^2)", ylabel="SNR", zlabel="Detection Range (km)", additional_params=additional_params)

    # 6. Multiple Scenarios Validation (Parametric wrt Instrument Aperture and Target Temperature)
    scenarios = []
    for temp, aperture in zip(TEMPERATURES, aperture_sizes):
        # Use realistic physics-based calculations to evaluate detection range
        spectral_radiance = BlackbodyRadiation.spectralRadiantExitance(temp, WAVELENGTH_BAND[0])
        effective_radiance = spectral_radiance * GrayBodyRadiation.calculate_emissivity(spectral_radiance, spectral_radiance)
        detection_range = DetectionPhysics.calculateDetectionRange(effective_radiance, sensorSensitivity=6e-19) * (aperture / 0.1)
        scenarios.append((time_values, detection_range, f'Temperature = {temp} K, Aperture = {aperture} m'))

    # Plot multiple scenarios to compare sensor performance with parametric variation
    PlotLib.plot_multiple_scenarios(scenarios, xlabel='Time (minutes)', ylabel='Detection Range (km)')

if __name__ == "__main__":
    main()
