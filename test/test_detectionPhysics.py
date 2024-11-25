from libraries.physicslib import DetectionPhysics

# Medium-wave IR sensor example
wavelength_mw_ir = 4e-6  # 4 µm in meters
aperture_diameter = 0.5  # meters

angular_resolution_mw_ir = DetectionPhysics.diffraction_limited_angular_resolution(wavelength_mw_ir, aperture_diameter)
print(f"Medium-Wave IR Sensor Angular Resolution: {angular_resolution_mw_ir * 1e6:.2f} µrad")  # Convert to µrad

power_emitted = 950  # W
nefd = 6.0e-19 *1e4     # W/cm²

snr_required = 10

detection_range = DetectionPhysics.detection_range(power_emitted, nefd, snr_required)
print(f"Detection Range: {detection_range / 1e3:.1f} km")  # Convert to km


area = 1.0  # m²
emissivity_visible = 0.4

reflectivity_area = DetectionPhysics.reflectivity_area_product(area, emissivity_visible)
print(f"Reflectivity-Area Product: {reflectivity_area:.2f} m²")

# Given parameters
diameter_sensor = 0.5  # m
range_to_target = 2.0e6  # m (2000 km)
power_emitted_target = 950  # W

power_received = DetectionPhysics.power_received_by_sensor(diameter_sensor, range_to_target, power_emitted_target)
print(f"Power Received by Sensor: {power_received:.2e} W")  # Expected 3.7e-12 W

temperature_target = 300  # K

emissivity_area = DetectionPhysics.emissivity_area_product(range_to_target, power_received, diameter_sensor, temperature_target)
print(f"Emissivity-Area Product: {emissivity_area:.2f} m²")  # Expected ~3.6 m²
