from poliastro.bodies import Earth
from astropy import units as u
from astropy.time import Time
from libraries.orbitlib import OrbitLib
from libraries.plotlib import OrbitPlotLib
import matplotlib.pyplot as plt

def plot_orbit_and_missile():
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
    missile_launch_altitude = Earth.R.to(u.km).value + 0.1  # Starting at 0.1 km (close to Earth's surface)
    missile_apogee = Earth.R.to(u.km).value + 200  # Apogee at 200 km altitude

    missile_orbit = OrbitLib.from_classical(
        name="Missile",
        attractor=Earth,
        a=(missile_launch_altitude + missile_apogee) / 2,  # Semi-major axis
        ecc=(missile_apogee - missile_launch_altitude) / (missile_apogee + missile_launch_altitude),  # Eccentricity
        inc=20.0,  # Inclination at 20 degrees
        raan=60.0,  # RAAN
        argp=0.0,  # Argument of perigee
        nu=0.0,  # True anomaly
        epoch=epoch
    )

    # Step 3: Plot the Orbits
    plt.figure(figsize=(10, 8))
    OrbitPlotLib.plot_orbit_3d(satellite_orbit.initial_orbit, name="Satellite Orbit")
    OrbitPlotLib.plot_orbit_3d(missile_orbit.initial_orbit, name="Missile Orbit")

    # Step 4: Configure and Show the Plot
    plt.title("Satellite and Missile Orbits")
    plt.xlabel("X Position (km)")
    plt.ylabel("Y Position (km)")
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    plot_orbit_and_missile()
