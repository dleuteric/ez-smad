# detectorlib.py

import numpy as np


class Detector:
    def __init__(self, quantum_efficiency, dark_current_params, shot_noise_params, thermal_noise_params, wavelength_min,
                 wavelength_max, f_number, netd):
        self.quantum_efficiency = quantum_efficiency
        self.dark_current_params = dark_current_params
        self.shot_noise_params = shot_noise_params
        self.thermal_noise_params = thermal_noise_params
        self.wavelength_min = wavelength_min  # Minimum wavelength for the detector band (in meters)
        self.wavelength_max = wavelength_max  # Maximum wavelength for the detector band (in meters)
        self.f_number = f_number  # Optical system f-number
        self.netd = netd  # Noise Equivalent Temperature Difference

    def calculate_dark_current(self, temperature):
        I0 = self.dark_current_params['I0']
        Eg = self.dark_current_params['Eg']
        k = self.dark_current_params['k']
        # Ensure the result is numeric
        return I0 * np.exp(-Eg / (k * temperature))

    def calculate_shot_noise(self, dark_current, bandwidth):
        q = self.shot_noise_params['q']
        # Validate that dark_current is a valid number
        if not isinstance(dark_current, (float, int)):
            raise ValueError(f"Expected dark_current to be a float, got {type(dark_current)}")

        return np.sqrt(2 * q * dark_current * bandwidth)

    def calculate_thermal_noise(self, temperature, bandwidth):
        k = self.thermal_noise_params['k']
        resistance = self.thermal_noise_params['resistance']
        return np.sqrt(4 * k * temperature * resistance * bandwidth)

    def calculate_nefd(self, aperture_size, integration_time, temperature, background_flux):
        # Calculate the bandwidth as the inverse of integration time
        bandwidth = 1 / integration_time

        # Calculate dark current
        dark_current = self.calculate_dark_current(temperature)

        # Calculate shot noise
        shot_noise = self.calculate_shot_noise(dark_current, bandwidth)

        # Calculate thermal noise
        thermal_noise = self.calculate_thermal_noise(temperature, bandwidth)

        # Total noise equivalent power (NEP)
        nep = np.sqrt(dark_current ** 2 + shot_noise ** 2 + thermal_noise ** 2)

        # Calculate NEFD (Noise Equivalent Flux Density)
        detector_area = np.pi * (aperture_size / 2) ** 2
        nefd = nep / (detector_area * bandwidth)

        return nefd


class SWIRDetector(Detector):
    def __init__(self):
        quantum_efficiency = 0.7
        dark_current_params = {'I0': 1e-9, 'Eg': 0.23, 'k': 1.38e-23}
        shot_noise_params = {'q': 1.6e-19}
        thermal_noise_params = {'k': 1.38e-23, 'resistance': 50}
        wavelength_min = 1e-6  # 1 µm
        wavelength_max = 3e-6  # 3 µm
        f_number = 2.5  # Placeholder value
        netd = 10e-3  # Placeholder value
        super().__init__(quantum_efficiency, dark_current_params, shot_noise_params, thermal_noise_params,
                         wavelength_min, wavelength_max, f_number, netd)


class MWIRDetector(Detector):
    def __init__(self):
        quantum_efficiency = 0.8
        dark_current_params = {'I0': 5e-9, 'Eg': 0.1, 'k': 1.38e-23}
        shot_noise_params = {'q': 1.6e-19}
        thermal_noise_params = {'k': 1.38e-23, 'resistance': 60}
        wavelength_min = 3e-6  # 3 µm
        wavelength_max = 8e-6  # 8 µm
        f_number = 2.0  # Placeholder value
        netd = 20e-3  # Placeholder value
        super().__init__(quantum_efficiency, dark_current_params, shot_noise_params, thermal_noise_params,
                         wavelength_min, wavelength_max, f_number, netd)


class LWIRDetector(Detector):
    def __init__(self):
        quantum_efficiency = 0.2
        dark_current_params = {'I0': 1e-6, 'Eg': 0.05, 'k': 1.38e-23}
        shot_noise_params = {'q': 1.6e-19}
        thermal_noise_params = {'k': 1.38e-23, 'resistance': 100}
        wavelength_min = 7.6e-6  # 7.6 µm
        wavelength_max = 9e-6  # 9 µm
        f_number = 2.0  # Value from datasheet
        netd = 0.05  # Value from datasheet
        super().__init__(quantum_efficiency, dark_current_params, shot_noise_params, thermal_noise_params,
                         wavelength_min, wavelength_max, f_number, netd)


class DetectorFactory:
    @staticmethod
    def create_detector(detector_type):
        if detector_type == "SWIR":
            return SWIRDetector()
        elif detector_type == "MWIR":
            return MWIRDetector()
        elif detector_type == "LWIR":
            return LWIRDetector()
        else:
            raise ValueError(f"Unknown detector type: {detector_type}")