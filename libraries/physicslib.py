# physicslib.py

import math

import numpy as np
from scipy import integrate
from astropy import units as u

# Physical Constants
PI = 3.141592653589793
PLANCK_CONSTANT = 6.62607015e-34        # Planck constant [J·s]
SPEED_OF_LIGHT = 299792458              # Speed of light in vacuum [m/s]
BOLTZMANN_CONSTANT = 1.380649e-23       # Boltzmann constant [J/K]
STEFAN_BOLTZMANN_CONSTANT = 5.670374419e-8  # Stefan-Boltzmann constant [W/m^2·K^4]
WIEN_DISPLACEMENT_CONSTANT = 2.897771955e-3 # Wien's displacement constant [m·K]

class BlackbodyRadiation:
    """
    Class for calculating blackbody radiation properties.

    Equations:
    - Wien's Displacement Law (λ_max [m] = b / T [K])
    - Stefan-Boltzmann Law (M [W/m^2] = σ * T^4)
    - Planck’s Law (Spectral Radiant Exitance M_λ [W/m^2·m^-1])
    """

    @staticmethod
    def spectral_radiant_exitance(T: float, wavelength_m: float) -> float:
        """
        Calculate spectral radiant exitance M_λ [W/m^2·m^-1] using Planck’s Law.

        M_λ = (2 * PI * h * c^2) / (λ^5 * (exp(h * c / (λ * k_B * T)) - 1))

        :param T: Temperature in Kelvin [K]
        :param wavelength_m: Wavelength in meters [m]
        :return: Spectral radiant exitance in W/m^2·m^-1
        """
        exponent = (PLANCK_CONSTANT * SPEED_OF_LIGHT) / (wavelength_m * BOLTZMANN_CONSTANT * T)
        numerator = 2 * PI * PLANCK_CONSTANT * SPEED_OF_LIGHT ** 2
        denominator = wavelength_m ** 5 * (math.exp(exponent) - 1)
        return numerator / denominator

    @staticmethod
    def peak_wavelength(T: float) -> float:
        """
        Calculate peak wavelength λ_max [m] using Wien's Displacement Law.

        λ_max = b / T

        :param T: Temperature in Kelvin [K]
        :return: Peak wavelength in meters [m]
        """
        return WIEN_DISPLACEMENT_CONSTANT / T

    @staticmethod
    def radiant_exitance(T: float) -> float:
        """
        Calculate total radiant exitance M [W/m^2] using Stefan-Boltzmann Law.

        M = σ * T^4

        :param T: Temperature in Kelvin [K]
        :return: Radiant exitance in watts per square meter [W/m^2]
        """
        return STEFAN_BOLTZMANN_CONSTANT * T ** 4


class GrayBodyRadiation:
    """
    Class for calculating gray body radiation properties.
    """

    @staticmethod
    def emissivity(object_radiance: float, blackbody_radiance: float) -> float:
        """
        Calculate emissivity ε(λ, T).

        ε(λ, T) = M_object(λ, T) / M_blackbody(λ, T)

        :param object_radiance: Radiance of the object [W/m^2·sr·m]
        :param blackbody_radiance: Radiance of a blackbody at the same temperature [W/m^2·sr·m]
        :return: Emissivity (dimensionless)
        """
        return object_radiance / blackbody_radiance

    @staticmethod
    def radiant_exitance(T: float, emissivity: float) -> float:
        """
        Calculate total radiant exitance M [W/m^2] for a gray body using Stefan-Boltzmann Law.

        M = ε * σ * T^4

        :param T: Temperature in Kelvin [K]
        :param emissivity: Emissivity of the gray body (dimensionless)
        :return: Radiant exitance in watts per square meter [W/m^2]
        """
        return emissivity * STEFAN_BOLTZMANN_CONSTANT * T ** 4

    @staticmethod
    def spectral_radiant_exitance(T: float, wavelength_m: float, emissivity: float) -> float:
        """
        Calculate spectral radiant exitance M_λ [W/m^2·m^-1] for a gray body.

        M_λ = ε * M_blackbody(λ, T)

        :param T: Temperature in Kelvin [K]
        :param wavelength_m: Wavelength in meters [m]
        :param emissivity: Emissivity of the gray body (dimensionless)
        :return: Spectral radiant exitance in W/m^2·m^-1
        """
        blackbody_exitance = BlackbodyRadiation.spectral_radiant_exitance(T, wavelength_m)
        return emissivity * blackbody_exitance

    @staticmethod
    def absorptivity(emissivity: float) -> float:
        """
        Calculate absorptivity α for an object in thermal equilibrium.

        α = ε

        :param emissivity: Emissivity of the gray body (dimensionless)
        :return: Absorptivity (dimensionless)
        """
        return emissivity

    @staticmethod
    def reflectivity(emissivity: float) -> float:
        """
        Calculate reflectivity ρ for an opaque object.

        ρ = 1 - ε

        :param emissivity: Emissivity of the gray body (dimensionless)
        :return: Reflectivity (dimensionless)
        """
        return 1 - emissivity

class EquilibriumTemperature:
    """
    Class for calculating the equilibrium temperature of objects in space.
    """

    @staticmethod
    def equilibrium_temperature(
        alpha_v: float,
        epsilon_ir: float,
        area_ratio: float = 0.25,
        solar_flux: float = 1360.0,
        albedo: float = 0.3,
        earth_ir_flux: float = 240.0,
        sigma: float = STEFAN_BOLTZMANN_CONSTANT
    ) -> float:
        """
        Calculate the equilibrium temperature T_eq [K] for an object in space.

        T_eq = [ (A_C / A_S) * ((S + S_R) * α_V + ε_IR * E) / (ε_IR * σ) ]^{1/4}

        :param alpha_v: Absorptivity averaged over visible and near-infrared band (dimensionless)
        :param epsilon_ir: Emissivity averaged over the infrared band (dimensionless)
        :param area_ratio: Ratio of cross-sectional area to surface area (A_C / A_S)
        :param solar_flux: Solar flux S [W/m^2]
        :param albedo: Earth's albedo (dimensionless)
        :param earth_ir_flux: Earth's infrared flux E [W/m^2]
        :param sigma: Stefan-Boltzmann constant [W/m^2·K^4]
        :return: Equilibrium temperature in Kelvin [K]
        """
        # Calculate reflected solar flux S_R
        s_r = albedo * solar_flux

        # Numerator: (S + S_R) * α_V + ε_IR * E
        numerator = ((solar_flux + s_r) * alpha_v) + (epsilon_ir * earth_ir_flux)

        # Denominator: ε_IR * σ
        denominator = epsilon_ir * sigma

        # Temperature calculation
        temp_fourth_power = (area_ratio * numerator) / denominator
        equilibrium_temp = temp_fourth_power ** 0.25

        return equilibrium_temp

    @staticmethod
    def equilibrium_temperature_in_shadow(
            epsilon_ir: float,
            earth_ir_flux: float = 240.0,
            sigma: float = STEFAN_BOLTZMANN_CONSTANT
    ) -> float:
        """
        Calculate the equilibrium temperature T_eq [K] for an object in Earth's shadow.

        T_eq = [ (A_C / A_S) * (ε_IR * E) / (ε_IR * σ) ]^{1/4}
             = [ (A_C / A_S) * E / σ ]^{1/4}

        Since (A_C / A_S) = 1/4 for a sphere, and ε_IR cancels out.

        :param epsilon_ir: Emissivity averaged over the infrared band (dimensionless)
        :param earth_ir_flux: Earth's infrared flux E [W/m^2]
        :param sigma: Stefan-Boltzmann constant [W/m^2·K^4]
        :return: Equilibrium temperature in Kelvin [K]
        """
        area_ratio = 0.25  # For a sphere
        numerator = area_ratio * earth_ir_flux
        denominator = sigma
        temp_fourth_power = numerator / denominator
        equilibrium_temp = temp_fourth_power ** 0.25
        return equilibrium_temp


class TransientThermalModel:
    """
    Class for modeling transient thermal behavior of objects in space.
    """

    @staticmethod
    def temperature_change(
        T: float,
        t: float,
        params: dict
    ) -> float:
        """
        Calculate the rate of change of temperature dT/dt at time t.

        :param T: Current temperature [K]
        :param t: Time [s] (not used in this model but included for ODE solvers)
        :param params: Dictionary containing required parameters
        :return: Rate of change of temperature dT/dt [K/s]
        """
        # Unpack parameters
        C = params['heat_capacity']          # Heat capacity [J/K]
        A_S = params['surface_area']         # Surface area [m^2]
        A_C = params['cross_sectional_area'] # Cross-sectional area [m^2]
        alpha_V = params['alpha_v']          # Visible absorptivity
        epsilon_IR = params['epsilon_ir']    # Infrared emissivity
        S = params['solar_flux']             # Solar flux [W/m^2]
        albedo = params['albedo']            # Earth's albedo
        E = params['earth_ir_flux']          # Earth's infrared flux [W/m^2]
        sigma = STEFAN_BOLTZMANN_CONSTANT    # Stefan-Boltzmann constant [W/m^2·K^4]

        # Calculate reflected solar flux
        S_R = albedo * S

        # Power absorbed
        P_A = A_C * ((S + S_R) * alpha_V + epsilon_IR * E)  # [W]

        # Power emitted
        P_E = A_S * epsilon_IR * sigma * T ** 4  # [W]

        # Rate of change of temperature
        dT_dt = (P_A - P_E) / C  # [K/s]

        return dT_dt


class DetectionPhysics:
    """
    Class for detection-related calculations, including sensor resolution and detection ranges.
    """

    @staticmethod
    def diffraction_limited_angular_resolution(wavelength_m: float, aperture_diameter_m: float) -> float:
        """
        Calculate the diffraction-limited angular resolution θ [radians].

        θ = λ / d

        :param wavelength_m: Wavelength of the radiation [m]
        :param aperture_diameter_m: Diameter of the sensor aperture [m]
        :return: Diffraction-limited angular resolution [radians]
        """
        return wavelength_m / aperture_diameter_m

    @staticmethod
    def detection_range_old_0(power_w: float, nefd_w_m2: float, snr: float) -> float:
        """
        Calculate the detection range R [m] for an object emitting power P [W].

        R = sqrt(P / (4π * NEFD * S/N))

        :param power_w: Power emitted by the object [W]
        :param nefd_w_m2: Noise Equivalent Flux Density [W/m^2]
        :param snr: Signal-to-noise ratio (dimensionless)
        :return: Detection range in meters [m]
        """
        return math.sqrt(power_w / (4 * PI * nefd_w_m2 * snr))


    @staticmethod
    def detection_range_old_1(
        power_emitted: float,
        nefd_w_m2: float,
        snr: float,
        band_radiant_exitance: float,
        quantum_efficiency: float = 0.7,
        optical_efficiency: float = 0.6
    ) -> float:
        """
        Calculate the detection range R [m] for an object emitting power P [W].

        R = sqrt[(P * QE * OE) / (4π * NEFD * S/N)]

        :param power_emitted: Power emitted by the object [W]
        :param nefd_w_m2: Noise Equivalent Flux Density [W/m^2]
        :param snr: Signal-to-noise ratio (dimensionless)
        :param band_radiant_exitance: Band-specific radiant exitance [W/m^2] (output from BandRadiation.band_radiant_exitance())
        :param quantum_efficiency: Efficiency of the detector in converting photons to signals (0.7 as default)
        :param optical_efficiency: Efficiency of the telescope optics to transmit radiation (default to 0.6)
        :return: Detection range in meters [m]
        """
        # Effective power contribution due to efficiencies and exitance
        effective_power_emitted = power_emitted * quantum_efficiency * optical_efficiency

        # Incorporate radiant exitance into the detection range calculation
        if band_radiant_exitance == 0:
            raise ValueError("Band radiant exitance must be greater than 0.")

        # Effective radiant power detected in band range
        effective_radiant_power = effective_power_emitted * band_radiant_exitance

        # Calculate detection range
        detection_range = math.sqrt(effective_radiant_power / (4 * PI * nefd_w_m2 * snr))

        return detection_range

    @staticmethod
    def detection_range(power_emitted, nefd_w_m2, snr, band_radiant_exitance, f_number, netd):
        """
        Calculates detection range based on the power emitted, NEFD, and signal-to-noise ratio (SNR).
        Incorporates additional parameters like f-number and NETD for more accuracy.

        Parameters:
        - power_emitted: Total emitted power (W)
        - nefd_w_m2: Noise Equivalent Flux Density (W/m²)
        - snr: Signal-to-noise ratio
        - band_radiant_exitance: Radiant exitance for the specific wavelength band (W/m²)
        - f_number: Cold shield f/# for the detector (affects light collection efficiency)
        - netd: Noise Equivalent Temperature Difference for the LWIR detector

        Returns:
        - Detection range (m)
        """
        # Validation checks
        if None in [power_emitted, nefd_w_m2, snr, band_radiant_exitance, f_number, netd]:
            raise ValueError("One or more parameters are None. Please check inputs.")
        # Existing calculations
        # Factor in the cold shield f-number in optical efficiency
        optical_efficiency_factor = 1 / (f_number ** 2)

        # Adjust the power received at the detector considering f-number effect
        power_received = power_emitted * optical_efficiency_factor

        # Calculate the noise power based on NEFD and NETD
        noise_power = nefd_w_m2 * netd

        # Use the SNR to determine if the detection range is adequate
        if noise_power == 0:
            return 0  # To prevent division by zero
        detection_range = (power_received / (nefd_w_m2 * snr)) ** 0.5

        return detection_range

    @staticmethod
    def reflectivity_area_product(area_m2: float, emissivity_v: float) -> float:
        """
        Calculate the reflectivity-area product ρV * A for a target in the visible band.

        ρV * A = (1 - εV) * A

        :param area_m2: Cross-sectional area of the target [m²]
        :param emissivity_v: Emissivity in the visible band (dimensionless)
        :return: Reflectivity-area product [m²]
        """
        return (1 - emissivity_v) * area_m2

    @staticmethod
    def power_received_by_sensor(diameter_m: float, range_m: float, power_emitted_w: float) -> float:
        """
        Calculate the power received by the sensor P_S [W].

        P_S = (π d^2 / 4) * (1 / R^2) * P_T

        :param diameter_m: Diameter of the sensor optics [m]
        :param range_m: Range to the target [m]
        :param power_emitted_w: Power emitted by the target [W]
        :return: Power received by the sensor [W]
        """
        return (PI * diameter_m ** 2 / 4) * (1 / range_m ** 2) * power_emitted_w

    @staticmethod
    def emissivity_area_product(range_m: float, power_received_w: float, diameter_m: float, temperature_k: float) -> float:
        """
        Calculate the emissivity-area product ε_IR A [m^2].

        ε_IR A = (4 R^2 P_S) / (σ π d^2 T^4)

        :param range_m: Range to the target [m]
        :param power_received_w: Power received by the sensor [W]
        :param diameter_m: Diameter of the sensor optics [m]
        :param temperature_k: Temperature of the target [K]
        :return: Emissivity-area product [m^2]
        """
        numerator = 4 * range_m**2 * power_received_w
        denominator = STEFAN_BOLTZMANN_CONSTANT * PI * diameter_m**2 * temperature_k**4
        return numerator / denominator


class BandRadiation:
    """
    Class for calculating band-specific radiant exitance using integration.
    """

    @staticmethod
    def band_radiant_exitance(T: float, wavelength_min: float, wavelength_max: float) -> float:
        """
        Calculate the band-specific radiant exitance for a graybody object by integrating Planck's law over a given wavelength range.

        M_band = ∫ [ (2 * PI * h * c^2) / (λ^5 * (exp(h * c / (λ * k_B * T)) - 1)) ] dλ

        :param T: Temperature in Kelvin [K]
        :param wavelength_min: Minimum wavelength of the band in meters [m]
        :param wavelength_max: Maximum wavelength of the band in meters [m]
        :return: Radiant exitance in W/m^2
        """

        # Reuse planck_spectral_radiance from BlackbodyRadiation
        def spectral_radiance(wavelength):
            return BlackbodyRadiation.spectral_radiant_exitance(T, wavelength)

        # Use numerical integration (quad) to integrate over the given range
        radiant_exitance, _ = integrate.quad(spectral_radiance, wavelength_min, wavelength_max)
        return radiant_exitance

class PhysicsLib:
    @staticmethod
    def thrust_acceleration(T, m):
        """
        Compute the thrust acceleration vector.

        :param T: Thrust vector [Tx, Ty, Tz] in Newtons
        :param m: Mass of the object in kg
        :return: Thrust acceleration vector in m/s^2
        """
        return np.array(T) / m

    @staticmethod
    def mass_rate_change(T, I_sp):
        """
        Compute the rate of change of mass due to thrust.

        :param T: Thrust vector [Tx, Ty, Tz] in Newtons
        :param I_sp: Specific impulse in seconds
        :return: Mass rate of change in kg/s
        """
        g0 = 9.80665  # Standard gravity in m/s^2
        T_magnitude = np.linalg.norm(T)
        return -T_magnitude / (g0 * I_sp)

    @staticmethod
    def atmospheric_density(h):
        """
        Compute atmospheric density as a function of altitude.
        """
        if h < 0:
            raise ValueError("Altitude cannot be negative.")
        rho_0 = 1.226  # Sea level density in kg/m^3
        scale_height = 7254.24  # Scale height in meters
        rho = rho_0 * np.exp(-h / scale_height)
        print(f"Altitude: {h}, Atmospheric Density: {rho}")
        return rho if rho > 1e-10 else 0  # Set a threshold to avoid extremely small values

    @staticmethod
    def drag_acceleration(v, h, C_D, S, m):
        """
        Compute drag acceleration.
        """
        if h < 0:
            raise ValueError("Altitude cannot be negative.")
        rho = PhysicsLib.atmospheric_density(h)
        v_mag = np.linalg.norm(v)
        v_mag = max(v_mag, 1e-3)  # Ensure velocity magnitude is non-zero
        drag = -0.5 * C_D * S * rho * v_mag * v / m
        drag_magnitude_limit = 10 * PhysicsLib.gravity_acceleration(np.array([0, 0, 6371000]))[2]  # Example limit
        drag = np.clip(drag, -drag_magnitude_limit, drag_magnitude_limit)  # Limit drag to realistic values
        print(f"Drag parameters -> rho: {rho}, v_mag: {v_mag}, drag: {drag}")
        return drag

    @staticmethod
    def gravity_acceleration(r):
        """
        Compute gravitational acceleration.

        :param r: Position vector in meters
        :return: Gravity acceleration vector in m/s^2
        """
        mu = 3.986004418e14  # Earth's gravitational parameter in m^3/s^2
        r_mag = np.linalg.norm(r)
        if r_mag < 1e-6:  # Add threshold to avoid numerical instability
            r_mag = 1e-6
        return -mu * r / r_mag ** 3

    @staticmethod
    def non_inertial_acceleration(v, r):
        """
        Compute non-inertial accelerations (Coriolis and centrifugal).

        :param v: Velocity vector in m/s
        :param r: Position vector in meters
        :return: Non-inertial acceleration vector in m/s^2
        """
        omega_e = np.array([0, 0, 7.2921150e-5])  # Earth's rotation rate in rad/s
        coriolis = -2 * np.cross(omega_e, v)
        centrifugal = -np.cross(omega_e, np.cross(omega_e, r))
        return coriolis + centrifugal