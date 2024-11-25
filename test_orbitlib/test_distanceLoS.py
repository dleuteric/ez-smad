from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from astropy import units as u
from astropy.time import Time
from libraries.orbitlib import OrbitLib
from libraries.mathlib import VectorOperations
import numpy as np

def main():
    # Define the epoch for both the satellite and missile
    epoch = Time("2025-01-01 12:00:00")

    # Step 1: Define the Satellite Orbit at 1000 km Altitude
    satellite_orbit = OrbitLib.from_classical(
        name="Satellite",
        attractor=Earth,
        a=Earth.R.to(u.km).value + 1000,  # Semi-major axis (Earth radius + 1000 km altitude)
        ecc=0.0,  # Circular orbit
        inc=45.0,  # Inclination at 45 degrees
        raan=30.0,  # RAAN
        argp=0.0,  # Argument of perigee
        nu=0.0,  # True anomaly
        epoch=epoch
    )

    # Step 2: Define the Missile Trajectory
    # Assuming ballistic missile-like trajectory (simplified as an orbit)
    missile_launch_altitude = Earth.R.to(u.km).value + 0.1  # Starting at 0.1 km (close to Earth's surface)
    missile_apogee = Earth.R.to(u.km).value + 200  # Apogee at 200 km altitude

    missile_orbit = OrbitLib.from_classical(
        name="Missile",
        attractor=Earth,
        a=(missile_launch_altitude + missile_apogee) / 2,  # Semi-major axis (average altitude between launch and apogee)
        ecc=(missile_apogee - missile_launch_altitude) / (missile_apogee + missile_launch_altitude),  # Approximate eccentricity
        inc=20.0,  # Inclination at 20 degrees
        raan=60.0,  # RAAN
        argp=0.0,  # Argument of perigee
        nu=0.0,  # True anomaly
        epoch=epoch
    )

    # Step 3: Propagate Both Orbits to the Same Time (e.g., after 10 minutes)
    time_delta = 600  # Propagate by 600 seconds (10 minutes)
    propagated_satellite = satellite_orbit.propagate(time_delta)
    propagated_missile = missile_orbit.propagate(time_delta)

    # Step 4: Retrieve Position Vectors
    satellite_position = propagated_satellite.r  # Position vector of the satellite
    missile_position = propagated_missile.r  # Position vector of the missile

    # Step 5: Compute the Distance Between Satellite and Missile
    distance = VectorOperations.distance(satellite_position, missile_position)

    # Step 6: Display the Results
    print(f"Satellite Position (after 10 minutes): {satellite_position} km")
    print(f"Missile Position (after 10 minutes): {missile_position} km")
    print(f"Distance between Satellite and Missile: {distance:.2f} km")

if __name__ == "__main__":
    main()
