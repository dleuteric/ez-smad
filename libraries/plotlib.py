import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from poliastro.plotting.static import StaticOrbitPlotter

# Set global defaults for color-blind friendly plots
plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 16,
    "legend.fontsize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "axes.prop_cycle": plt.cycler(color=["#0072B2", "#E69F00", "#56B4E9", "#F0E442", "#009E73", "#D55E00", "#CC79A7"])
})

class PlotLib:
    """
    Library for centralizing all plot-related functions for visualizing the different quantities in the SBIRS-Low analysis.
    """

    @staticmethod
    def set_plot_defaults():
        """
        Set the default styles for the plots.
        """
        plt.style.use('ggplot')  # Use a default style to make plots visually appealing.

    @staticmethod
    def plot_relationship(x, y, xlabel, ylabel, title, legend=None, log_scale=False):
        """
        Generalized function to plot the relationship between any two parameters.
        :param x: Array-like, data for x-axis.
        :param y: Array-like, data for y-axis.
        :param xlabel: String, label for x-axis.
        :param ylabel: String, label for y-axis.
        :param title: String, title of the plot.
        :param legend: String, optional legend for the data.
        :param log_scale: Boolean, whether to use log scale for the y-axis.
        """
        plt.figure(figsize=(10, 6))
        plt.plot(x, y, marker='o', linestyle='-', linewidth=2)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        if log_scale:
            plt.yscale('log')
        if legend:
            plt.legend([legend])
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    @staticmethod
    def plot_detection_range_vs_snr_parametric(snr_values, detection_ranges, temperatures, fixed_parameter='cross section', parameter_value=None):
        """
        Plot detection range as a function of SNR for different object temperatures or cross sections.
        :param snr_values: Array-like, SNR values.
        :param detection_ranges: 2D list or array, detection ranges for different temperatures or cross sections.
        :param temperatures: List, temperatures or cross sections for different series in the plot.
        :param fixed_parameter: String, indicating the fixed parameter ('temperature' or 'cross section').
        :param parameter_value: The value of the fixed parameter.
        """
        plt.figure(figsize=(10, 6))
        for idx, temp in enumerate(temperatures):
            plt.plot(snr_values, detection_ranges[idx], marker='o', linestyle='-', linewidth=2, label=f'{fixed_parameter.capitalize()} = {temp}')
        plt.xlabel('Signal-to-Noise Ratio (SNR)')
        plt.ylabel('Detection Range (km)')
        plt.title('Detection Range vs SNR for Different Conditions')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        note = f"Fixed Aperture Size = 0.3 m, Fixed {fixed_parameter.capitalize()} = {parameter_value}"
        plt.figtext(0.5, -0.1, note, wrap=True, horizontalalignment='center', fontsize=10)
        plt.show()

    @staticmethod
    def plot_detection_range_vs_aperture_parametric(aperture_sizes, detection_ranges, fixed_param_value, fixed_param_name='SNR', fixed_temp=230):
        """
        Plot detection range as a function of aperture size, with fixed values of SNR and temperature.
        :param aperture_sizes: Array-like, aperture sizes of the sensor.
        :param detection_ranges: Array-like, detection ranges for each aperture size.
        :param fixed_param_value: Value of the fixed parameter (e.g., SNR or temperature).
        :param fixed_param_name: String, specifying which parameter is fixed.
        :param fixed_temp: Temperature at which data is evaluated.
        """
        plt.figure(figsize=(10, 6))
        plt.plot(aperture_sizes, detection_ranges, marker='o', linestyle='-', linewidth=2)
        plt.xlabel('Aperture Size (m)')
        plt.ylabel('Detection Range (km)')
        plt.title('Detection Range vs Aperture Size for Fixed Conditions')
        plt.grid(True)
        plt.tight_layout()
        note = f"Fixed Temperature = {fixed_temp} K, Fixed {fixed_param_name} = {fixed_param_value}"
        plt.figtext(0.5, -0.1, note, wrap=True, horizontalalignment='center', fontsize=10)
        plt.show()

    @staticmethod
    def plot_time_evolution_of_temperature(time_values, temperature_values, material_name, total_time_minutes=120):
        """
        Plot the temperature evolution over time for different materials.
        :param time_values: Array-like, time points.
        :param temperature_values: Array-like, temperature values over time.
        :param material_name: String, name of the material for annotation.
        :param total_time_minutes: Total time in minutes for the temperature evolution.
        """
        plt.figure(figsize=(10, 6))
        plt.plot(time_values, temperature_values, marker='o', linestyle='-', linewidth=2, label=f'Material: {material_name}')
        plt.xlabel('Time (minutes)')
        plt.ylabel('Temperature (K)')
        plt.title('Temperature Evolution Over Time')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        note = f"Assuming a total simulation time of {total_time_minutes} minutes. This duration may reflect orbital geometry considerations."
        plt.figtext(0.5, -0.1, note, wrap=True, horizontalalignment='center', fontsize=10)
        plt.show()

    @staticmethod
    def plot_parametric_analysis(param1, param2, z_values, xlabel, ylabel, zlabel, additional_params=None):
        """
        Create a 3D surface plot to show the trade-off between two parameters with additional parameters for more insight.
        :param param1: Array-like, first independent parameter.
        :param param2: Array-like, second independent parameter.
        :param z_values: 2D array, dependent variable.
        :param xlabel: String, label for x-axis.
        :param ylabel: String, label for y-axis.
        :param zlabel: String, label for z-axis.
        :param additional_params: Dictionary of additional parameters and their values.
        """
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Use indexing='ij' to ensure the shapes match
        param1_mesh, param2_mesh = np.meshgrid(param1, param2, indexing='ij')

        # Now param1_mesh, param2_mesh, and z_values all have the same shape
        ax.plot_surface(param1_mesh, param2_mesh, z_values, cmap='viridis', alpha=0.7)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_zlabel(zlabel)
        ax.set_title('Parametric Analysis')
        if additional_params:
            note = ', '.join([f"{key} = {value}" for key, value in additional_params.items()])
            plt.figtext(0.5, -0.05, note, wrap=True, horizontalalignment='center', fontsize=10)
        plt.tight_layout()
        plt.show()

    @staticmethod
    def plot_multiple_scenarios(scenarios, xlabel, ylabel):
        """
        Plot multiple scenario curves on the same figure.
        :param scenarios: List of tuples, each containing x and y data for different scenarios.
        :param xlabel: String, label for x-axis.
        :param ylabel: String, label for y-axis.
        """
        plt.figure(figsize=(10, 6))
        for scenario in scenarios:
            x_data, y_data, label = scenario
            plt.plot(x_data, y_data, marker='o', linestyle='-', linewidth=2, label=label)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title('Comparison of Multiple Scenarios')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    @staticmethod
    def plot_multiple_relationships(x, y_dict, xlabel, ylabel, title):
        """
        Plot multiple relationships in a single plot for comparison.
        :param x: Array-like, data for x-axis.
        :param y_dict: Dictionary of y-data with labels as keys.
        :param xlabel: String, label for x-axis.
        :param ylabel: String, label for y-axis.
        :param title: String, title of the plot.
        """
        plt.figure(figsize=(10, 6))
        for label, y_data in y_dict.items():
            plt.plot(x, y_data, marker='o', linestyle='-', linewidth=2, label=label)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    @staticmethod
    def plot_band_radiant_exitance(temperatures, results):
        """
        Plot radiant exitance as a function of temperature for different wavelength bands.
        :param temperatures: List of temperatures in Kelvin.
        :param results: Dictionary of radiant exitance values, with temperature as keys and wavebands as subkeys.
        """
        plt.figure(figsize=(10, 6))

        # Prepare data for plotting
        bands = list(next(iter(results.values())).keys())
        for band in bands:
            y_values = [results[temp][band] for temp in temperatures]
            plt.plot(temperatures, y_values, marker='o', linestyle='-', linewidth=2, label=band)

        # Plotting details
        plt.xlabel('Temperature (K)')
        plt.ylabel('Radiant Exitance (W/m^2)')
        plt.title('Band-Specific Radiant Exitance vs Temperature')
        plt.legend(title="Wavebands")
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    @staticmethod
    def set_plot_defaults():
        plt.style.use("seaborn")
        plt.rcParams.update({
            "axes.grid": True,
            "grid.alpha": 0.6,
            "figure.figsize": (10, 6),
            "lines.linewidth": 2,
            "lines.markersize": 8,
        })

    @staticmethod
    def plot_multiple_detection_ranges_vs_temperature(range_results, xlabel="Temperature (K)",
                                                      ylabel="Detection Range (km)",
                                                      title="Detection Range vs Temperature"):
        """
        Plot detection range vs temperature for different detector configurations.

        :param range_results: A dictionary where keys are labels for each combination, and values are lists of detection ranges.
        :param xlabel: Label for the x-axis.
        :param ylabel: Label for the y-axis.
        :param title: Title for the plot.
        """
        plt.figure()
        for label, detection_ranges in range_results.items():
            temperatures = range_results[label]['temperatures']
            detection_range = range_results[label]['detection_ranges']
            plt.plot(temperatures, detection_range, marker='o', linestyle='-', label=label)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    @staticmethod
    def plot_comparison_between_wavebands(detection_data, xlabel="Temperature (K)", ylabel="Detection Range (km)",
                                          title="Detection Range Comparison Across Wavebands"):
        """
        Plot comparison of detection ranges across different wavebands.

        :param detection_data: A dictionary where keys are wavebands and values are dicts containing temperature, detection range.
        :param xlabel: Label for the x-axis.
        :param ylabel: Label for the y-axis.
        :param title: Title for the plot.
        """
        plt.figure()
        for band, data in detection_data.items():
            plt.plot(data["temperatures"], data["detection_ranges"], marker='o', linestyle='-', label=f"{band}")
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    @staticmethod
    def plot_detection_range_comparison(temperatures, range_results, aperture_sizes, snr, xlabel="Temperature (K)", ylabel="Detection Range (km)", title="Detection Range vs Temperature for Different Aperture Sizes"):
        """
        Plot detection ranges for different aperture sizes, for a given SNR.
        :param temperatures: Array-like, temperature values.
        :param range_results: Dictionary containing detection range results for different aperture sizes.
        :param aperture_sizes: List of aperture sizes.
        :param snr: The signal-to-noise ratio for this plot.
        :param xlabel: String, label for x-axis.
        :param ylabel: String, label for y-axis.
        :param title: String, title of the plot.
        """
        plt.figure(figsize=(12, 8))
        for aperture, detection_range in range_results.items():
            plt.plot(temperatures, detection_range, marker='o', linestyle='-', linewidth=2, label=f"Aperture Size: {aperture} m")
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(f"{title} (SNR = {snr})")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

class OrbitPlotLib:
    @staticmethod
    def plot_orbit_2d(orbit, name):
        """
        Plot the orbit in a 2D plane using matplotlib.
        """
        fig, ax = plt.subplots(figsize=(10, 10))
        plotter = StaticOrbitPlotter(ax)
        plotter.plot(orbit, label=name)
        plt.title(f"Orbit of {name}")
        plt.xlabel("x (km)")
        plt.ylabel("y (km)")
        plt.grid(True)
        plt.show()

    @staticmethod
    def plot_orbit_3d(orbit, name):
        """
        Plot the orbit in a 3D plane using matplotlib.
        """
        from mpl_toolkits.mplot3d import Axes3D  # Import here to avoid unnecessary dependency if 3D is not used
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')
        sample_points = orbit.sample(100)  # Get sample points along the orbit
        ax.plot(sample_points.x.value, sample_points.y.value, sample_points.z.value, label=name)
        ax.set_title(f"3D Orbit of {name}")
        ax.set_xlabel("x (km)")
        ax.set_ylabel("y (km)")
        ax.set_zlabel("z (km)")
        ax.legend()
        plt.show()

    @staticmethod
    def plot_multiple_orbits_2d(orbits, labels):
        """
        Plot multiple orbits in a 2D plane using matplotlib.
        :param orbits: List of Orbit objects to plot.
        :param labels: List of labels corresponding to each orbit.
        """
        fig, ax = plt.subplots(figsize=(10, 10))
        plotter = StaticOrbitPlotter(ax)
        for orbit, label in zip(orbits, labels):
            plotter.plot(orbit, label=label)
        plt.title("Multiple Orbits in 2D")
        plt.xlabel("x (km)")
        plt.ylabel("y (km)")
        plt.grid(True)
        plt.legend()
        plt.show()

    @staticmethod
    def plot_ballistic_trajectory(solution, label="Ballistic Trajectory"):
        """
        Plot a 3D ballistic trajectory using the solution from solve_ivp.

        Parameters:
        - solution: Solution object from solve_ivp containing trajectory data
        - label: Label for the plot
        """
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')

        x = solution.y[0]
        y = solution.y[1]
        z = solution.y[2]

        ax.plot(x, y, z, label=label)
        ax.set_title("3D Ballistic Trajectory")
        ax.set_xlabel("X Position (km)")
        ax.set_ylabel("Y Position (km)")
        ax.set_zlabel("Z Position (km)")
        ax.legend()
        plt.show()
